@testset "Spatial Partitioning" begin
    @testset "Cell Indexing" begin
        # nx = 10, ny = 10, cell_size = 5.0 -> inv_cell_size = 0.2
        inv_cell_size = 0.2f0
        @test cell_index(0.1f0,  0.1f0,  inv_cell_size, 10, 10) == 1
        @test cell_index(5.1f0,  0.1f0,  inv_cell_size, 10, 10) == 2
        @test cell_index(0.1f0,  5.1f0,  inv_cell_size, 10, 10) == 11
        @test cell_index(49.9f0, 49.9f0, inv_cell_size, 10, 10) == 100

        # Out-of-bounds must be clamped (not throw)
        @test cell_index(-1.0f0, -1.0f0, inv_cell_size, 10, 10) == 1
        @test cell_index(60.0f0, 60.0f0, inv_cell_size, 10, 10) == 100
    end

    @testset "Neighbor Table" begin
        table = EcoEvolutionSim.build_neighbor_table(10, 10)
        @test length(table.count) == 100

        # Bottom-left corner (cell 1): 4 neighbors
        @test sort(table.neighbors[1:table.count[1], 1]) == [1, 2, 11, 12]

        # Top-right corner (cell 100): 4 neighbors
        @test sort(table.neighbors[1:table.count[100], 100]) == [89, 90, 99, 100]

        # Interior cell (cell 15 -> ix=5, iy=2): full 3x3 = 9 neighbors
        @test sort(table.neighbors[1:table.count[15], 15]) == [4, 5, 6, 14, 15, 16, 24, 25, 26]

        # Edge cells have fewer than 9 neighbors
        @test table.count[1]  == 4   # corner
        @test table.count[2]  == 6   # top edge (not corner)
        @test table.count[15] == 9   # interior
    end

    # -----------------------------------------------------------------------
    # Helper: build a tiny sim for grid tests
    # -----------------------------------------------------------------------
    function make_grid_sim(N, ncells, nx, ny, cell_size)
        traits = (test = zeros(Float32, N),)
        agents = Agents(N, N, Float32(nx * cell_size), Float32(ny * cell_size), traits)
        grid = CellGrid(ncells, N)
        nt = EcoEvolutionSim.build_neighbor_table(nx, ny)
        env = EnvironmentState(ncells, nx, ny, Float32(cell_size), 1f0 / cell_size, grid, nt)
        cfg = Config(world_width=Float32(nx * cell_size), world_length=Float32(ny * cell_size), cell_size=Float32(cell_size))
        return Simulation(cfg, agents, env, x -> nothing, ())
    end

    @testset "Grid Building — single thread" begin
        # world 10×10, cell_size 5 -> 2×2 grid, 4 cells
        # cell 1: x∈[0,5), y∈[0,5)  -> agents 1,2
        # cell 2: x∈[5,10),y∈[0,5)  -> agent 3
        # cell 3: x∈[0,5), y∈[5,10) -> (none)
        # cell 4: x∈[5,10),y∈[5,10) -> agent 4
        # agent 5 is dead
        sim = make_grid_sim(5, 4, 2, 2, 5)
        agents = sim.agents
        agents.x .= Float32[1.0, 4.0, 6.0, 8.0, 1.0]
        agents.y .= Float32[1.0, 4.0, 2.0, 7.0, 1.0]
        agents.alive .= [true, true, true, true, false]
        agents.n_alive = 4

        build_cell_grid!(sim)

        grid = sim.env.grid
        @test grid.cell_count == [2, 1, 0, 1]
        @test grid.cell_start == [1, 3, 4, 4]
        @test sum(grid.cell_count) == 4   # dead agent excluded

        c1 = sort(grid.agent_index[grid.cell_start[1] : grid.cell_start[1] + grid.cell_count[1] - 1])
        @test c1 == [1, 2]
        @test grid.agent_index[grid.cell_start[2]] == 3
        @test grid.agent_index[grid.cell_start[4]] == 4
    end

    @testset "Grid Building — idempotency (rebuild gives same result)" begin
        sim = make_grid_sim(10, 4, 2, 2, 5)
        sim.agents.x .= rand(Float32, 10) .* 10f0
        sim.agents.y .= rand(Float32, 10) .* 10f0
        build_cell_grid!(sim)
        counts1 = copy(sim.env.grid.cell_count)
        starts1 = copy(sim.env.grid.cell_start)
        build_cell_grid!(sim)
        @test sim.env.grid.cell_count == counts1
        @test sim.env.grid.cell_start == starts1
    end

    @testset "Grid Building — no living agents gives all-zero counts" begin
        sim = make_grid_sim(5, 4, 2, 2, 5)
        sim.agents.alive .= false
        sim.agents.n_alive = 0
        build_cell_grid!(sim)
        @test all(sim.env.grid.cell_count .== 0)
    end

    @testset "nearby_agents — distance filtering & self-exclusion" begin
        # Place agents at known positions in a 10×10 world, cell_size=5
        sim = make_grid_sim(5, 4, 2, 2, 5)
        agents = sim.agents
        # Agent 1 at (2,2); agents 2-4 at various distances; agent 5 dead
        agents.x .= Float32[2.0, 2.5, 5.0, 8.0, 2.0]
        agents.y .= Float32[2.0, 2.5, 2.0, 8.0, 2.0]
        agents.alive .= [true, true, true, true, false]
        agents.n_alive = 4
        build_cell_grid!(sim)

        # r=1.0: agent 2 is sqrt(0.5)≈0.71 away -> included; agent 3 is 3.0 away -> excluded
        neighbors = collect(nearby_agents(sim, 1, 1.0f0))
        @test 2 ∈ neighbors     # close
        @test 1 ∉ neighbors     # self excluded
        @test 3 ∉ neighbors     # too far
        @test 5 ∉ neighbors     # dead

        # r=5.0: agents 2 and 3 within range
        neighbors_wide = collect(nearby_agents(sim, 1, 5.0f0))
        @test 2 ∈ neighbors_wide
        @test 3 ∈ neighbors_wide
        @test 5 ∉ neighbors_wide  # dead stays excluded
    end

    @testset "nearby_agents — empty neighbourhood" begin
        sim = make_grid_sim(3, 4, 2, 2, 5)
        agents = sim.agents
        agents.x .= Float32[1.0, 9.0, 9.0]
        agents.y .= Float32[1.0, 9.0, 9.0]
        agents.alive .= [true, false, false]
        agents.n_alive = 1
        build_cell_grid!(sim)
        @test isempty(collect(nearby_agents(sim, 1, 2.0f0)))
    end
end
