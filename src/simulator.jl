using Base.Threads
using Random
using TOML

############################################################
# Configuration and Models
############################################################

@enum MovementStrategy begin
    RANDOM_WALK = 1
    LANGEVIN = 2
    CORRELATED_RW = 3
    ACTIVE_BROWNIAN = 4
end

struct Config
    # Simulation settings
    n_agents::Int
    world_size::Float32
    cell_size::Float32
    
    # Movement settings
    movement_strategy::MovementStrategy
    base_speed::Float32
    noise_strength::Float32
    friction::Float32 # For Langevin
    persistence::Float32 # For Correlated RW
    
    # Ecological settings
    interaction_radius2::Float32
    competition_sigma::Float32
    mating_sigma::Float32
    predation_threshold::Float32
    
    # Trait indices (placeholder for expansion)
    # 1: niche/resource preference
    # ...
end

function Config(;
    n_agents = 500,
    world_size = 50.0f0,
    cell_size = 5.0f0,
    movement_strategy = RANDOM_WALK,
    base_speed = 0.1f0,
    noise_strength = 0.1f0,
    friction = 0.1f0,
    persistence = 0.8f0,
    interaction_radius2 = 1.0f0,
    competition_sigma = 0.2f0,
    mating_sigma = 0.05f0,
    predation_threshold = 0.3f0
)
    return Config(
        n_agents, world_size, cell_size,
        movement_strategy, base_speed, noise_strength, friction, persistence,
        interaction_radius2, competition_sigma, mating_sigma, predation_threshold
    )
end

function load_config(path::String)
    data = TOML.parsefile(path)
    
    # Map string strategy to enum
    strategy_map = Dict(
        "RANDOM_WALK" => RANDOM_WALK,
        "LANGEVIN" => LANGEVIN,
        "CORRELATED_RW" => CORRELATED_RW,
        "ACTIVE_BROWNIAN" => ACTIVE_BROWNIAN
    )
    
    sim = data["simulation"]
    mov = data["movement"]
    eco = data["ecology"]
    
    return Config(
        n_agents = Int(sim["n_agents"]),
        world_size = Float32(sim["world_size"]),
        cell_size = Float32(sim["cell_size"]),
        movement_strategy = strategy_map[mov["strategy"]],
        base_speed = Float32(mov["base_speed"]),
        noise_strength = Float32(mov["noise_strength"]),
        friction = Float32(mov["friction"]),
        persistence = Float32(mov["persistence"]),
        interaction_radius2 = Float32(eco["interaction_radius2"]),
        competition_sigma = Float32(eco["competition_sigma"]),
        mating_sigma = Float32(eco["mating_sigma"]),
        predation_threshold = Float32(eco["predation_threshold"])
    )
end

############################################################
# Agent storage (Structure of Arrays)
############################################################

struct Agents
    x::Vector{Float32}
    y::Vector{Float32}
    
    # Movement state
    vx::Vector{Float32}
    vy::Vector{Float32}
    theta::Vector{Float32}
    
    # Ecological state
    trait::Vector{Float32}
    energy::Vector{Float32}
end

function Agents(N, world_size)
    return Agents(
        rand(Float32, N) .* world_size,
        rand(Float32, N) .* world_size,
        zeros(Float32, N),
        zeros(Float32, N),
        rand(Float32, N) .* 2f0 .* Float32(pi),
        rand(Float32, N),
        ones(Float32, N)
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
# Simulation object
############################################################

struct Simulation
    config::Config
    agents::Agents
    grid::CellGrid
    neighbor_cells::Vector{Vector{Int}}
    
    # Helper fields for convenience
    nx::Int
    ny::Int
    ncells::Int
end

############################################################
# Convert position → grid cell
############################################################

@inline function cell_index(x, y, cell_size, nx)
    ix = Int(floor(x / cell_size)) + 1
    iy = Int(floor(y / cell_size)) + 1
    return ix + nx * (iy - 1)
end

############################################################
# Precompute neighbor cell table
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

############################################################
# Build spatial cell list
############################################################

function build_cell_grid!(sim::Simulation)
    agents = sim.agents
    grid = sim.grid
    ncells = sim.ncells
    cell_size = sim.config.cell_size
    nx = sim.nx

    for c in 1:ncells
        grid.cell_count[c][] = 0
        grid.cell_offset[c][] = 0
    end

    N = length(agents.x)

    Threads.@threads for i in 1:N
        c = cell_index(agents.x[i], agents.y[i], cell_size, nx)
        Threads.atomic_add!(grid.cell_count[c], 1)
    end

    s = 1
    for c in 1:ncells
        grid.cell_start[c] = s
        s += grid.cell_count[c][]
    end

    Threads.@threads for i in 1:N
        c = cell_index(agents.x[i], agents.y[i], cell_size, nx)
        idx = Threads.atomic_add!(grid.cell_offset[c], 1)
        grid.agent_index[grid.cell_start[c] + idx] = i
    end
end

############################################################
# Ecological interaction kernel
############################################################

@inline function interaction_kernel!(i, j, sim)
    agents = sim.agents
    config = sim.config

    ti = agents.trait[i]
    tj = agents.trait[j]

    d = ti - tj
    ad = abs(d)

    # Competition
    comp = exp(-(ad^2) / (2f0 * config.competition_sigma^2))
    agents.energy[i] -= 0.01f0 * comp
    agents.energy[j] -= 0.01f0 * comp

    # Mating
    mate = exp(-(ad^2) / (2f0 * config.mating_sigma^2))
    if mate > 0.8f0
        agents.energy[i] += 0.05f0
        agents.energy[j] += 0.05f0
    end

    # Predation
    if d > config.predation_threshold
        agents.energy[i] += 0.1f0
        agents.energy[j] -= 0.1f0
    elseif d < -config.predation_threshold
        agents.energy[j] += 0.1f0
        agents.energy[i] -= 0.1f0
    end
end

############################################################
# Interaction step
############################################################

function interaction_step!(sim::Simulation)
    agents = sim.agents
    grid = sim.grid
    config = sim.config
    N = length(agents.x)

    Threads.@threads for i in 1:N
        xi = agents.x[i]
        yi = agents.y[i]
        c = cell_index(xi, yi, config.cell_size, sim.nx)

        for nc in sim.neighbor_cells[c]
            start = grid.cell_start[nc]
            stop = start + grid.cell_count[nc][] - 1
            for k in start:stop
                j = grid.agent_index[k]
                if j ≤ i continue end

                dx = agents.x[j] - xi
                dy = agents.y[j] - yi
                r2 = dx * dx + dy * dy

                if r2 < config.interaction_radius2
                    interaction_kernel!(i, j, sim)
                end
            end
        end
    end
end

############################################################
# Movement step
############################################################

function movement_step!(sim::Simulation)
    agents = sim.agents
    config = sim.config
    world = config.world_size
    N = length(agents.x)

    Threads.@threads for i in 1:N
        if config.movement_strategy == RANDOM_WALK
            agents.x[i] += config.base_speed * (rand(Float32) - 0.5f0)
            agents.y[i] += config.base_speed * (rand(Float32) - 0.5f0)

        elseif config.movement_strategy == LANGEVIN
            # dv = -friction * v * dt + sqrt(2*D*dt) * dW
            # Here base_speed acts as noise strength for simplicity
            agents.vx[i] = agents.vx[i] * (1f0 - config.friction) + config.noise_strength * (rand(Float32) - 0.5f0)
            agents.vy[i] = agents.vy[i] * (1f0 - config.friction) + config.noise_strength * (rand(Float32) - 0.5f0)
            agents.x[i] += agents.vx[i]
            agents.y[i] += agents.vy[i]

        elseif config.movement_strategy == CORRELATED_RW
            # Change angle slightly
            agents.theta[i] += config.noise_strength * (rand(Float32) - 0.5f0)
            agents.x[i] += config.base_speed * cos(agents.theta[i])
            agents.y[i] += config.base_speed * sin(agents.theta[i])

        elseif config.movement_strategy == ACTIVE_BROWNIAN
            # Constant speed, rotational diffusion
            agents.theta[i] += config.noise_strength * (rand(Float32) - 0.5f0)
            agents.x[i] += config.base_speed * cos(agents.theta[i])
            agents.y[i] += config.base_speed * sin(agents.theta[i])
            # (In this simple version AB is similar to CRW but usually speed is fixed)
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
    interaction_step!(sim)
    movement_step!(sim)
end

############################################################
# Initialization
############################################################

function init_simulation(config::Config)
    nx = Int(config.world_size / config.cell_size)
    ny = nx
    ncells = nx * ny

    agents = Agents(config.n_agents, config.world_size)
    neighbor_table = build_neighbor_table(nx, ny)
    grid = CellGrid(ncells, config.n_agents)

    return Simulation(
        config,
        agents,
        grid,
        neighbor_table,
        nx,
        ny,
        ncells
    )
end

# Backward compatibility or convenience
function init_simulation(N, world_size, cell_size)
    config = Config(n_agents=N, world_size=world_size, cell_size=cell_size)
    return init_simulation(config)
end
