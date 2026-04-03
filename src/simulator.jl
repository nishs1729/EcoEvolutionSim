# Initialize simulation
function init_simulation(config_path::String = "config.toml", traits_path::String = "script/traits.toml")
    @timeit to "initialization" begin
        config = load_config(config_path)
        nx = Int(config.world_size / config.cell_size)
        ny = nx
        ncells = nx * ny

        trait_specs = load_traits(traits_path)
        traits = initialize_traits(trait_specs, config.n_agents)
        agents = Agents(config.n_agents, config.world_size, traits)
        neighbor_table = build_neighbor_table(nx, ny)
        grid = CellGrid(ncells, config.n_agents)
        inv_cell_size = 1.0f0 / config.cell_size
        env = EnvironmentState(ncells, nx, ny, config.cell_size, inv_cell_size, grid, neighbor_table)
        movement_kernel = select_movement_kernel(config.strategy)

        return Simulation(config, agents, env, movement_kernel)
    end
end

# One simulation step
function step!(sim::Simulation)
    @timeit to "step" begin
        @timeit to "build_grid" build_cell_grid!(sim)
        ecology_step!(sim)
        @timeit to "movement" movement_step!(sim)
    end
end
