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

@option struct SimConfig
    n_agents::Int = 500
    world_size::Float32 = 50.0f0
    cell_size::Float32 = 5.0f0
end

@option struct MovConfig
    strategy::MovementStrategy = RANDOM_WALK
    base_speed::Float32 = 0.1f0
    noise_strength::Float32 = 0.1f0
    friction::Float32 = 0.1f0 # For Langevin
    persistence::Float32 = 0.8f0 # For Correlated RW
end

@option struct EcoConfig
    interaction_radius2::Float32 = 1.0f0
    competition_sigma::Float32 = 0.2f0
    mating_sigma::Float32 = 0.05f0
    predation_threshold::Float32 = 0.3f0
end

@option struct Config
    sim::SimConfig = SimConfig()
    mov::MovConfig = MovConfig()
    eco::EcoConfig = EcoConfig()
    traits::Dict{String, Dict{String,Any}} = Dict{String, Any}("resource_preference" => Dict("type" => "continuous"))
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

############################################################
# Agent storage (SoA)
############################################################

struct Agents
    x::Vector{Float32}
    y::Vector{Float32}
    vx::Vector{Float32}
    vy::Vector{Float32}
    theta::Vector{Float32}
    energy::Vector{Float32}
    traits::Dict{String, Vector{Float32}}
end

function Agents(N, world_size, traits_dict::Dict{String, Vector{Float32}})
    # For backward compatibility, pick the first trait or use a default
    default_trait = if haskey(traits_dict, "resource_preference")
        traits_dict["resource_preference"]
    elseif !isempty(traits_dict)
        first(values(traits_dict))
    else
        rand(Float32, N)
    end

    return Agents(
        rand(Float32, N) .* world_size,
        rand(Float32, N) .* world_size,
        zeros(Float32, N),
        zeros(Float32, N),
        rand(Float32, N) .* 2f0 .* Float32(pi),
        ones(Float32, N),
        traits_dict,
        default_trait
    )
end

############################################################
# Spatial grid
############################################################

struct CellGrid
    cell_start::Vector{Int}
    cell_count::Vector{Threads.Atomic{Int}}
    cell_offset::Vector{Threads.Atomic{Int}}
    agent_index::Vector{Int}
end

function CellGrid(ncells, N)
    return CellGrid(
        zeros(Int, ncells),
        [Threads.Atomic{Int}(0) for _ in 1:ncells],
        [Threads.Atomic{Int}(0) for _ in 1:ncells],
        Vector{Int}(undef, N)
    )
end

############################################################
# Environment state
############################################################

include("environment.jl")
using .Environment

############################################################
# Registries
############################################################

using .FitnessRegistry
using .TraitRegistry

############################################################
# Simulation object
############################################################

struct Simulation
    config::Config
    agents::Agents
    env::EnvironmentState
    fitness_fn::Function
end

include("fitness.jl")
using .Fitness
include("traits.jl")
using .Traits

############################################################
# Utility: Position → grid cell
############################################################

@inline function cell_index(x, y, cell_size, nx)
    ix = Int(floor(x / cell_size)) + 1
    iy = Int(floor(y / cell_size)) + 1
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
        c = cell_index(agents.x[i], agents.y[i], cell_size, nx)
        Threads.atomic_add!(env.grid.cell_count[c], 1)
    end

    s = 1
    for c in 1:ncells
        env.grid.cell_start[c] = s
        s += env.grid.cell_count[c][]
    end

    Threads.@threads for i in 1:N
        c = cell_index(agents.x[i], agents.y[i], cell_size, nx)
        idx = Threads.atomic_add!(env.grid.cell_offset[c], 1)
        env.grid.agent_index[env.grid.cell_start[c] + idx] = i
    end
end

############################################################
# Movement step
############################################################

function movement_step!(sim::Simulation)
    agents = sim.agents
    config = sim.config
    world = config.sim.world_size
    N = length(agents.x)

    Threads.@threads for i in 1:N
        if config.mov.strategy == RANDOM_WALK
            agents.x[i] += config.mov.base_speed * (rand(Float32) - 0.5f0)
            agents.y[i] += config.mov.base_speed * (rand(Float32) - 0.5f0)

        elseif config.mov.strategy == LANGEVIN
            agents.vx[i] = agents.vx[i] * (1f0 - config.mov.friction) + config.mov.noise_strength * (rand(Float32) - 0.5f0)
            agents.vy[i] = agents.vy[i] * (1f0 - config.mov.friction) + config.mov.noise_strength * (rand(Float32) - 0.5f0)
            agents.x[i] += agents.vx[i]
            agents.y[i] += agents.vy[i]

        elseif config.mov.strategy == CORRELATED_RW
            agents.theta[i] += config.mov.noise_strength * (rand(Float32) - 0.5f0)
            agents.x[i] += config.mov.base_speed * cos(agents.theta[i])
            agents.y[i] += config.mov.base_speed * sin(agents.theta[i])

        elseif config.mov.strategy == ACTIVE_BROWNIAN
            agents.theta[i] += config.mov.noise_strength * (rand(Float32) - 0.5f0)
            agents.x[i] += config.mov.base_speed * cos(agents.theta[i])
            agents.y[i] += config.mov.base_speed * sin(agents.theta[i])
        end

        # Reflective boundaries
        if agents.x[i] < 0.0f0
            agents.x[i] = -agents.x[i]
            agents.vx[i] = -agents.vx[i]
            agents.theta[i] = Float32(pi) - agents.theta[i]
        elseif agents.x[i] > world
            agents.x[i] = 2.0f0 * world - agents.x[i]
            agents.vx[i] = -agents.vx[i]
            agents.theta[i] = Float32(pi) - agents.theta[i]
        end

        if agents.y[i] < 0.0f0
            agents.y[i] = -agents.y[i]
            agents.vy[i] = -agents.vy[i]
            agents.theta[i] = -agents.theta[i]
        elseif agents.y[i] > world
            agents.y[i] = 2.0f0 * world - agents.y[i]
            agents.vy[i] = -agents.vy[i]
            agents.theta[i] = -agents.theta[i]
        end

        agents.x[i] = clamp(agents.x[i], 0.0f0, world)
        agents.y[i] = clamp(agents.y[i], 0.0f0, world)
    end
end

############################################################
# One simulation step
############################################################

function step!(sim::Simulation)
    build_cell_grid!(sim)
    compute_interactions!(sim, sim.fitness_fn)
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

function init_simulation(config::Config)
    nx = Int(config.sim.world_size / config.sim.cell_size)
    ny = nx
    ncells = nx * ny

    # Initialize traits from registry
    traits_dict = Dict{String, Vector{Float32}}()
    for (trait_name, trait_config) in config.traits
        trait_type = trait_config["type"]
        trait_constructor = TRAIT_REGISTRY[trait_type]
        # For 'bounded', pass extra args if they exist
        if trait_type == "bounded"
            min_val = Float32(get(trait_config, "min", 0.0f0))
            max_val = Float32(get(trait_config, "max", 1.0f0))
            traits_dict[trait_name] = trait_constructor(config.sim.n_agents, min_val, max_val)
        else
            traits_dict[trait_name] = trait_constructor(config.sim.n_agents)
        end
    end

    agents = Agents(config.sim.n_agents, config.sim.world_size, traits_dict)
    neighbor_table = build_neighbor_table(nx, ny)
    grid = CellGrid(ncells, config.sim.n_agents)
    
    env = EnvironmentState(config, grid, neighbor_table, nx, ny, ncells)

    fitness_name = "default_interaction"
    fitness_constructor = FITNESS_REGISTRY[fitness_name]
    fitness_fn = fitness_constructor(config)

    return Simulation(config, agents, env, fitness_fn)
end

function init_simulation(N, world_size, cell_size)
    config = Config(sim=SimConfig(n_agents=N, world_size=world_size, cell_size=cell_size))
    return init_simulation(config)
end
