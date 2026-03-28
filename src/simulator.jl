using Base.Threads
using Random
using TOML
using Configurations

include("config.jl")
include("agents.jl")
include("traits.jl")
include("space.jl")
include("movement.jl")

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

# One simulation step
function step!(sim::Simulation)
    build_cell_grid!(sim)
    # compute_interactions!(sim, sim.fitness_fn)
    age_step!(sim)
    death_step!(sim)
    movement_step!(sim)
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
