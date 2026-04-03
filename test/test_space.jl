@testset "Spatial Partitioning" begin
    @testset "Cell Indexing" begin
        # nx = 10, ny = 10, cell_size = 5.0 -> inv_cell_size = 0.2
        # cells are 1 to 100
        # 0.0 < x <= 5.0 -> ix = 1
        inv_cell_size = 0.2f0
        @test cell_index(0.1f0, 0.1f0, inv_cell_size, 10, 10) == 1
        @test cell_index(5.1f0, 0.1f0, inv_cell_size, 10, 10) == 2
        @test cell_index(0.1f0, 5.1f0, inv_cell_size, 10, 10) == 11
        @test cell_index(49.9f0, 49.9f0, inv_cell_size, 10, 10) == 100
        
        # out of bounds should be clamped
        @test cell_index(-1.0f0, -1.0f0, inv_cell_size, 10, 10) == 1
        @test cell_index(60.0f0, 60.0f0, inv_cell_size, 10, 10) == 100
    end
    
    @testset "Neighbor Table" begin
        table = EcoEvolutionSim.build_neighbor_table(10, 10)
        @test length(table.count) == 100
        
        # Bottom-left corner (cell 1)
        # Neighbors: (1,1), (2,1), (1,2), (2,2) -> cells 1, 2, 11, 12
        @test sort(table.neighbors[1:table.count[1], 1]) == [1, 2, 11, 12]
        
        # Top-right corner (cell 100)
        # Neighbors: (9,9), (10,9), (9,10), (10,10) -> cells 89, 90, 99, 100
        @test sort(table.neighbors[1:table.count[100], 100]) == [89, 90, 99, 100]
        
        # Middle cell (cell 15 -> x=5, y=2)
        # Neighbors: (4,1),(5,1),(6,1), (4,2),(5,2),(6,2), (4,3),(5,3),(6,3)
        # -> cells 4,5,6, 14,15,16, 24,25,26
        @test sort(table.neighbors[1:table.count[15], 15]) == [4, 5, 6, 14, 15, 16, 24, 25, 26]
    end
    
    @testset "Grid Building" begin
        N = 5
        traits = (test = zeros(Float32, N),)
        agents = Agents(N, N, 10.0f0, 10.0f0, traits)
        
        # Fix coordinates for predictable cells
        # world 10.0, cell_size 5.0 -> nx=2, ny=2, ncells=4
        # cell 1: [0,5)x[0,5) -> agents 1, 2
        # cell 2: [5,10)x[0,5) -> agent 3
        # cell 3: [0,5)x[5,10) -> none
        # cell 4: [5,10)x[5,10) -> agent 4
        # agent 5 is dead
        agents.x .= Float32[1.0, 4.0, 6.0, 8.0, 1.0]
        agents.y .= Float32[1.0, 4.0, 2.0, 7.0, 1.0]
        agents.alive .= [true, true, true, true, false]
        
        grid = CellGrid(4, N)
        neighbor_table = NeighborTable(Matrix{Int32}(undef, 0, 0), Int8[])
        env = EnvironmentState(4, 2, 2, 5.0f0, 0.2f0, grid, neighbor_table)
        sim = Simulation(Config(world_width=10.0f0, world_length=10.0f0, cell_size=5.0f0), agents, env, x -> nothing)
        
        build_cell_grid!(sim)
        
        @test grid.cell_count == [2, 1, 0, 1]
        @test grid.cell_start == [1, 3, 4, 4]
        
        # Verify agents are in correct cells
        # cell 1 start at 1, length 2
        c1_agents = sort(grid.agent_index[grid.cell_start[1]:grid.cell_start[1]+grid.cell_count[1]-1])
        @test c1_agents == [1, 2]
        
        # cell 2 start at 3, length 1
        c2_agents = grid.agent_index[grid.cell_start[2]:grid.cell_start[2]+grid.cell_count[2]-1]
        @test c2_agents == [3]
        
        # cell 4 start at 4, length 1
        c4_agents = grid.agent_index[grid.cell_start[4]:grid.cell_start[4]+grid.cell_count[4]-1]
        @test c4_agents == [4]
        
        # Agent 5 is dead and should not be counted
        @test sum(grid.cell_count) == 4
    end
end
