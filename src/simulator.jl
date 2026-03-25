using Base.Threads
using Random
using TOML
using Configurations

############################################################
# Configuration and Models
############################################################

@enum MovementStrategy begin
    RANDOM_WALK = 1
    LANGEVIN = 2
    CORRELATED_RW = 3
    ACTIVE_BROWNIAN = 4
end

@option struct Config
    # Simulation
    n_agents::Int = 500
    world_size::Float32 = 50.0f0
    cell_size::Float32 = 5.0f0
    steps_per_frame::Int = 10

    # Movement
    strategy::MovementStrategy = RANDOM_WALK
    base_speed::Float32 = 0.1f0

    # Ecology / interactions
    r_interact::Float32 = 1.0f0
end

function Configurations.from_dict(::Type{MovementStrategy}, x::String)
    s = Symbol(uppercase(x))
    for i in instances(MovementStrategy)
        if Symbol(i) == s
            return i
        end
    end
    error("'$x' is not a valid MovementStrategy. Options: $(instances(MovementStrategy))")
end

Base.convert(::Type{MovementStrategy}, x::String) = Configurations.from_dict(MovementStrategy, x)

function load_config(path::String)
    return from_toml(Config, path)
end

struct TraitSpec
    mean::Float32
    sigma::Float32
end

TRAIT_SPECS = Dict(
    "dispersal" => TraitSpec(1.0f0, 0.2f0),
    "fecundity" => TraitSpec(5.0f0, 1.0f0)
)

############################################################
# Agent storage (SoA)
############################################################

struct Agents
    x::Vector{Float32}
    y::Vector{Float32}
    vx::Vector{Float32}
    vy::Vector{Float32}
    gender::Vector{UInt8} # 0: female, 1: male
    energy::Vector{Float32}
    traits::Dict{String, Vector{Float32}}
end

function Agents(N, world_size, traits_dict::Dict{String, Vector{Float32}})
    return Agents(
        rand(Float32, N) .* world_size,  # x
        rand(Float32, N) .* world_size,  # y
        zeros(Float32, N),  # vx
        zeros(Float32, N),  # vy
        rand(UInt8(0):UInt8(1), N),  # gender
        rand(Float32, N),  # energy
        traits_dict  # traits
    )
end

############################################################
# Spatial grid
############################################################

struct CellGrid
    cell_start::Vector{Int}
    cell_count::Vector{Threads.Atomic{Int}}
    cell_offset::Vector{Threads.Atomic{Int}}
    agent_index::Vector{Int32}
end

function CellGrid(ncells, N)
    return CellGrid(
        zeros(Int, ncells),
        [Threads.Atomic{Int}(0) for _ in 1:ncells],
        [Threads.Atomic{Int}(0) for _ in 1:ncells],
        Vector{Int32}(undef, N)
    )
end

############################################################
# Simulation object
############################################################

struct EnvironmentState
    ncells::Int
    nx::Int
    ny::Int
    cell_size::Float32
    grid::CellGrid
    neighbors::Vector{Vector{Int}}
end

struct Simulation{F}
    config::Config
    agents::Agents
    env::EnvironmentState
    movement_kernel!::F
end

############################################################
# Utility: Position → grid cell
############################################################

@inline function cell_index(x, y, cell_size, nx, ny)
    ix = clamp(Int(floor(x / cell_size)) + 1, 1, nx)
    iy = clamp(Int(floor(y / cell_size)) + 1, 1, ny)
    return ix + nx * (iy - 1)
end

############################################################
# Build spatial cell list
############################################################

function build_cell_grid!(sim::Simulation)
    agents = sim.agents
    env = sim.env
    ncells = env.ncells
    cell_size = env.cell_size
    nx = env.nx

    for c in 1:ncells
        env.grid.cell_count[c][] = 0
        env.grid.cell_offset[c][] = 0
    end

    N = length(agents.x)

    Threads.@threads for i in 1:N
        c = cell_index(agents.x[i], agents.y[i], cell_size, nx, env.ny)
        Threads.atomic_add!(env.grid.cell_count[c], 1)
    end

    s = 1
    for c in 1:ncells
        env.grid.cell_start[c] = s
        s += env.grid.cell_count[c][]
    end

    Threads.@threads for i in 1:N
        c = cell_index(agents.x[i], agents.y[i], cell_size, nx, env.ny)
        idx = Threads.atomic_add!(env.grid.cell_offset[c], 1)
        env.grid.agent_index[env.grid.cell_start[c] + idx] = i
    end
end

############################################################
# Movement step
############################################################
function select_movement_kernel(strategy::MovementStrategy)
    if strategy == RANDOM_WALK
        return movement_random_walk!
    elseif strategy == LANGEVIN
        return movement_langevin!
    elseif strategy == CORRELATED_RW
        return movement_correlated_rw!
    elseif strategy == ACTIVE_BROWNIAN
        return movement_active_brownian!
    else
        error("Unknown movement strategy: $strategy")
    end
end

function movement_step!(sim::Simulation)
    sim.movement_kernel!(sim)
end

@inline function reflect!(x, v, world)
    if x < 0f0
        return -x, -v
    elseif x > world
        return 2f0*world - x, -v
    end
    return x, v
end

function movement_random_walk!(sim::Simulation)

    agents = sim.agents
    config = sim.config
    world = config.world_size
    speed = config.base_speed
    half = 0.5f0 * speed

    x = agents.x
    y = agents.y
    vx = agents.vx
    vy = agents.vy

    Threads.@threads for i in eachindex(x)
        rng = Random.default_rng()
        @inbounds begin
            xi = x[i]
            yi = y[i]
            vxi = vx[i]
            vyi = vy[i]

            xi += rand(rng,Float32)*speed - half
            yi += rand(rng,Float32)*speed - half

            xi, vxi = reflect!(xi, vxi, world)
            yi, vyi = reflect!(yi, vyi, world)

            x[i] = xi
            y[i] = yi
            vx[i] = vxi
            vy[i] = vyi
        end
    end
end

function movement_langevin!(sim::Simulation)

    agents = sim.agents
    config = sim.config

    world = config.world_size
    noise = config.noise_strength
    friction = config.friction

    damp = 1f0 - friction
    halfnoise = 0.5f0 * noise

    x = agents.x
    y = agents.y
    vx = agents.vx
    vy = agents.vy

    Threads.@threads for i in eachindex(x)
        rng = Random.default_rng()
        @inbounds begin

            xi = x[i]
            yi = y[i]
            vxi = vx[i]
            vyi = vy[i]

            vxi = vxi*damp + noise*rand(rng,Float32) - halfnoise
            vyi = vyi*damp + noise*rand(rng,Float32) - halfnoise

            xi += vxi
            yi += vyi

            xi, vxi = reflect!(xi, vxi, world)
            yi, vyi = reflect!(yi, vyi, world)

            x[i] = xi
            y[i] = yi
            vx[i] = vxi
            vy[i] = vyi
        end
    end
end

function movement_correlated_rw!(sim::Simulation)

    agents = sim.agents
    config = sim.config

    world = config.world_size
    noise = config.noise_strength
    speed = config.base_speed

    halfnoise = 0.5f0 * noise
    π32 = Float32(π)
    twopi = 2f0 * π32

    x = agents.x
    y = agents.y
    theta = agents.theta

    Threads.@threads for i in eachindex(x)
        rng = Random.default_rng()
        @inbounds begin

            xi = x[i]
            yi = y[i]
            θ = theta[i]

            θ += noise*rand(rng,Float32) - halfnoise

            s,c = sincos(θ)
            xi += speed*c
            yi += speed*s

            if xi < 0f0 || xi > world
                θ = π32 - θ
            end

            if yi < 0f0 || yi > world
                θ = -θ
            end

            theta[i] = mod(θ,twopi)
            x[i] = xi
            y[i] = yi
        end
    end
end

function movement_active_brownian!(sim::Simulation)

    agents = sim.agents
    config = sim.config

    world = config.world_size
    noise = config.noise_strength
    speed = config.base_speed

    π32 = Float32(π)
    twopi = 2f0 * π32

    x = agents.x
    y = agents.y
    theta = agents.theta

    Threads.@threads for i in eachindex(x)
        rng = Random.default_rng()
        @inbounds begin

            xi = x[i]
            yi = y[i]
            θ = theta[i]

            θ += noise * randn(rng,Float32)

            s,c = sincos(θ)
            xi += speed*c
            yi += speed*s

            if xi < 0f0 || xi > world
                θ = π32 - θ
            end

            if yi < 0f0 || yi > world
                θ = -θ
            end

            theta[i] = mod(θ,twopi)
            x[i] = xi
            y[i] = yi
        end
    end
end

############################################################
# One simulation step
############################################################

function step!(sim::Simulation)
    build_cell_grid!(sim)
    # compute_interactions!(sim, sim.fitness_fn)
    movement_step!(sim)
end

############################################################
# Initialization
############################################################

function build_neighbor_table(nx, ny)
    ncells = nx * ny
    table = Vector{Vector{Int}}(undef, ncells)
    for c in 1:ncells
        iy = (c - 1) ÷ nx + 1
        ix = c - nx * (iy - 1)
        list = Int[]
        for dy in -1:1
            for dx in -1:1
                x2 = ix + dx
                y2 = iy + dy
                if 1 ≤ x2 ≤ nx && 1 ≤ y2 ≤ ny
                    push!(list, x2 + nx * (y2 - 1))
                end
            end
        end
        table[c] = list
    end
    return table
end

function initialize_traits(specs::Dict{String, TraitSpec}, n_agents)
    traits = Dict{String, Vector{Float32}}()
    buffer = Vector{Float32}(undef, n_agents)

    for (name, spec) in specs
        randn!(buffer)
        v = similar(buffer)
        @inbounds @simd for i in eachindex(buffer)
            v[i] = spec.mean + spec.sigma*buffer[i]
        end
        traits[name] = v
    end

    return traits
end

function init_simulation(config_path::String = "config.toml")
    config = load_config(config_path)
    nx = Int(config.world_size / config.cell_size)
    ny = nx
    ncells = nx * ny

    traits_dict = initialize_traits(TRAIT_SPECS, config.n_agents)
    agents = Agents(config.n_agents, config.world_size, traits_dict)
    neighbor_table = build_neighbor_table(nx, ny)
    grid = CellGrid(ncells, config.n_agents)
    env = EnvironmentState(ncells, nx, ny, config.cell_size, grid, neighbor_table)
    movement_kernel = select_movement_kernel(config.strategy)

    return Simulation(config, agents, env, movement_kernel)
end
