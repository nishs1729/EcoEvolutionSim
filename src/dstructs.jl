struct TraitSpec
    mean::Float32
    sigma::Float32
end

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

struct Agents{T}
    x::Vector{Float32}
    y::Vector{Float32}
    vx::Vector{Float32}
    vy::Vector{Float32}
    gender::Vector{UInt8} # 0: female, 1: male
    energy::Vector{Float32}
    age::Vector{Float32}
    alive::Vector{Bool}
    traits::T
end

struct CellGrid
    cell_start::Vector{Int}
    cell_count::Vector{Int}
    cell_offset::Vector{Int}
    agent_index::Vector{Int32}
end

struct EnvironmentState
    ncells::Int
    nx::Int
    ny::Int
    cell_size::Float32
    grid::CellGrid
    neighbors::Vector{Vector{Int}}
end

struct Simulation{F, T}
    config::Config
    agents::Agents{T}
    env::EnvironmentState
    movement_kernel!::F
end