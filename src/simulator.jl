# Initialize simulation
function init_simulation(config_path::String = "config.toml", traits_path::String = "script/traits.toml"; interactions::Tuple = ())
    @timeit to "initialization" begin
        config = load_config(config_path)
        nx = ceil(Int, config.world_width / config.cell_size)
        ny = ceil(Int, config.world_length / config.cell_size)
        ncells = nx * ny

        trait_specs = load_traits(traits_path)
        traits = initialize_traits(trait_specs, config.max_agents, config.n_agents)
        agents = Agents(config.max_agents, config.n_agents, config.world_width, config.world_length, traits)
        neighbor_table = build_neighbor_table(nx, ny)
        grid = CellGrid(ncells, config.max_agents)
        inv_cell_size = 1.0f0 / config.cell_size
        env = EnvironmentState(ncells, nx, ny, config.cell_size, inv_cell_size, grid, neighbor_table)
        movement_kernel = select_movement_kernel(config.strategy)

        return Simulation(config, agents, env, movement_kernel, interactions)
    end
end

# One simulation step
function step!(sim::Simulation)
    @timeit to "step" begin
        @timeit to "build_grid" build_cell_grid!(sim)
        interactions_step!(sim)
        ecology_step!(sim)
        @timeit to "movement" movement_step!(sim)
    end
end
