using Base.Threads
using Random
using TOML
using Configurations

include("config.jl")
include("agents.jl")
include("traits.jl")
include("space.jl")
include("movement.jl")

# Initialize simulation
function init_simulation(config_path::String = "config.toml", traits_path::String = "script/traits.toml")
    config = load_config(config_path)
    nx = Int(config.world_size / config.cell_size)
    ny = nx
    ncells = nx * ny
    
    trait_specs = load_traits(traits_path)
    traits = initialize_traits(trait_specs, config.n_agents)
    agents = Agents(config.n_agents, config.world_size, traits)
    neighbor_table = build_neighbor_table(nx, ny)
    grid = CellGrid(ncells, config.n_agents)
    env = EnvironmentState(ncells, nx, ny, config.cell_size, grid, neighbor_table)
    movement_kernel = select_movement_kernel(config.strategy)
    
    return Simulation(config, agents, env, movement_kernel)
end

# One simulation step
function step!(sim::Simulation)
    build_cell_grid!(sim)
    # compute_interactions!(sim, sim.fitness_fn)
    age_step!(sim)
    death_step!(sim)
    movement_step!(sim)
end