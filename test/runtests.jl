using EcoEvolutionSim
using Test
using Random

@testset "EcoEvolutionSim.jl" begin

    @testset "Initialization" begin
        N = 100
        world_size = 50.0f0
        cell_size = 5.0f0
        sim = init_simulation(N, world_size, cell_size)
        
        @test length(sim.agents.x) == N
        @test length(sim.agents.y) == N
        @test length(sim.agents.trait) == N
        @test length(sim.agents.energy) == N
        
        @test sim.world_size == world_size
        @test sim.cell_size == cell_size
        @test sim.nx == 10
        @test sim.ny == 10
        @test sim.ncells == 100
        
        @test all(0.0 .<= sim.agents.x .<= world_size)
        @test all(0.0 .<= sim.agents.y .<= world_size)
        @test all(0.0 .<= sim.agents.trait .<= 1.0)
        @test all(sim.agents.energy .== 1.0)
    end

    @testset "Spatial Grid Logic" begin
        # Test cell_index
        # Cell (1,1) is index 1
        @test cell_index(2.5f0, 2.5f0, 5.0f0, 10) == 1
        # Cell (2,1) is index 2
        @test cell_index(7.5f0, 2.5f0, 5.0f0, 10) == 2
        # Cell (1,2) is index 11 (nx=10)
        @test cell_index(2.5f0, 7.5f0, 5.0f0, 10) == 11
        
        # Test neighbor table for a center cell
        nx, ny = 10, 10
        table = EcoEvolutionSim.build_neighbor_table(nx, ny)
        # Center cell (5,5) -> index 45
        c = 5 + 10*(5-1) # 45
        neighbors = table[c]
        @test length(neighbors) == 9
        @test c in neighbors
        
        # Test neighbor table for a corner cell (1,1) -> index 1
        neighbors_corner = table[1]
        @test length(neighbors_corner) == 4 # (1,1), (2,1), (1,2), (2,2)
        @test 1 in neighbors_corner
        @test 2 in neighbors_corner
        @test 11 in neighbors_corner
        @test 12 in neighbors_corner
    end

    @testset "Movement and Boundaries" begin
        # Test reflective boundaries
        N = 1
        world_size = 10.0f0
        cell_size = 5.0f0
        sim = init_simulation(N, world_size, cell_size)
        
        # Manually set agent near boundary
        sim.agents.x[1] = 0.05f0
        sim.agents.y[1] = 0.05f0
        
        # We need to reach into internal functions or just run step!
        # Let's test movement_step! directly if possible
        # Since it's not exported, we can use EcoEvolutionSim.movement_step!
        
        # Mocking large movement to force reflection
        # movement_step! adds 0.1*(rand-0.5), so max step is 0.05
        # To test reflection, we can manually trigger it or adjust agent pos
        
        for _ in 1:100
            EcoEvolutionSim.movement_step!(sim)
            @test 0.0 <= sim.agents.x[1] <= world_size
            @test 0.0 <= sim.agents.y[1] <= world_size
        end
    end

    @testset "Interaction Mechanics" begin
        N = 2
        world_size = 10.0f0
        cell_size = 10.0f0 # Single cell
        sim = init_simulation(N, world_size, cell_size)
        
        # Place agents close together
        sim.agents.x[1], sim.agents.y[1] = 5.0f0, 5.0f0
        sim.agents.x[2], sim.agents.y[2] = 5.1f0, 5.1f0
        
        # Scenario 1: Identical traits (Competition)
        sim.agents.trait[1] = 0.5f0
        sim.agents.trait[2] = 0.5f0
        sim.agents.energy .= 1.0f0
        
        EcoEvolutionSim.build_cell_grid!(sim)
        EcoEvolutionSim.interaction_step!(sim)
        
        # energy should increase due to mating bonus (0.05) > competition (0.01)
        @test sim.agents.energy[1] > 1.0f0
        @test sim.agents.energy[2] > 1.0f0
        
        # Scenario 2: Predation
        # trait[1] = 0.8, trait[2] = 0.2
        # d = 0.6 > predation_threshold (0.3)
        sim.agents.trait[1] = 0.8f0
        sim.agents.trait[2] = 0.2f0
        sim.agents.energy .= 1.0f0
        
        EcoEvolutionSim.build_cell_grid!(sim)
        EcoEvolutionSim.interaction_step!(sim)
        
        # Agent 1 (predator) energy should increase
        # Agent 2 (prey) energy should decrease
        # Note: competition also applies, but predation is 0.1 and competition is ~0.01*exp(...)
        @test sim.agents.energy[1] > 1.0f0
        @test sim.agents.energy[2] < 1.0f0
        
        # Scenario 3: Mating
        # trait diff is small (0.01), mating width is 0.05
        sim.agents.trait[1] = 0.5f0
        sim.agents.trait[2] = 0.51f0
        sim.agents.energy .= 1.0f0
        
        EcoEvolutionSim.build_cell_grid!(sim)
        EcoEvolutionSim.interaction_step!(sim)
        
        # Mate bonus (0.05) > Competition cost (0.01)
        @test sim.agents.energy[1] > 1.0f0
        @test sim.agents.energy[2] > 1.0f0
    end

    @testset "Full Simulation Step" begin
        N = 100
        sim = init_simulation(N, 50.0f0, 5.0f0)
        
        # Check if step! runs without error
        try
            step!(sim)
            @test true
        catch e
            @test false
        end
    end

end
