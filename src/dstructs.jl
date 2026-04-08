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
    n_agents::Int = 500 # initial number
    max_agents::Int = 100_000 # max agents for memory allocation
    world_width::Float32 = 50.0f0 # km
    world_length::Float32 = 100.0f0 # km
    cell_size::Float32 = 5.0f0 # km
    steps_per_frame::Int = 10
    dt::Float32 = 1.0f0 # days
    
    # Movement
    strategy::MovementStrategy = RANDOM_WALK
    base_speed::Float32 = 0.1f0 # km/day
    noise_strength::Float32 = 0.1f0 # rad/sqrt(day) or km/day^(3/2)
    friction::Float32 = 0.1f0 # 1/day
    
    # Ecology / interactions
    r_interact::Float32 = 1.0f0 # km
    metabolic_rate::Float32 = 10.0f0 # Joules/day
end

mutable struct Agents{T}
    x::Vector{Float32}
    y::Vector{Float32}
    vx::Vector{Float32}
    vy::Vector{Float32}
    theta::Vector{Float32}           # heading angle (radians) — used by CORRELATED_RW & ACTIVE_BROWNIAN
    gender::Vector{UInt8}            # 0: female, 1: male
    energy::Vector{Float32}
    age::Vector{Float32}
    last_mating::Vector{Float32}
    alive::Vector{Bool}
    traits::T
    max_id::Int
    n_alive::Int                     # O(1) live population counter
    free_indices::Vector{Int}
    death_buffer::Vector{Vector{Int}}
end

struct NeighborTable
    neighbors::Matrix{Int32}
    count::Vector{Int8}
end

struct CellGrid
    cell_start::Vector{Int32}
    cell_count::Vector{Int32}
    cell_offset::Vector{Int32}
    agent_index::Vector{Int32}
    thread_counts::Matrix{Int32}   # ncells × max_tid — reused scratch, avoids per-step alloc
    thread_offsets::Matrix{Int32}  # ncells × max_tid — reused scratch, avoids per-step alloc
end

struct EnvironmentState
    ncells::Int
    nx::Int
    ny::Int
    cell_size::Float32
    inv_cell_size::Float32
    grid::CellGrid
    neighbors::NeighborTable
end

struct Simulation{F, T, I}
    config::Config
    agents::Agents{T}
    env::EnvironmentState
    movement_kernel!::F
    interactions::I
end