using EcoEvolutionSim
using Test
using Random

@testset "EcoEvolutionSim.jl" begin

    @testset "Configuration" begin
        config = Config(n_agents=200, world_size=100.0f0, movement_strategy=EcoEvolutionSim.LANGEVIN)
        @test config.n_agents == 200
        @test config.world_size == 100.0f0
        @test config.movement_strategy == EcoEvolutionSim.LANGEVIN
        
        # Test loading from file
        config_path = joinpath(@__DIR__, "..", "config.toml")
        loaded = load_config(config_path)
        @test loaded.n_agents == 500
        @test loaded.world_size == 50.0f0
        @test loaded.movement_strategy == EcoEvolutionSim.LANGEVIN
    end

    @testset "Initialization" begin
        N = 100
        world_size = 50.0f0
        cell_size = 5.0f0
        sim = init_simulation(N, world_size, cell_size)
        
        @test length(sim.agents.x) == N
        @test length(sim.agents.y) == N
        @test length(sim.agents.trait) == N
        @test length(sim.agents.energy) == N
        
        @test sim.config.world_size == world_size
        @test sim.config.cell_size == cell_size
        @test sim.env.nx == 10
        @test sim.env.ny == 10
        @test sim.env.ncells == 100
        
        @test all(0.0 .<= sim.agents.x .<= world_size)
        @test all(0.0 .<= sim.agents.y .<= world_size)
        @test all(0.0 .<= sim.agents.trait .<= 1.0)
        @test all(sim.agents.energy .== 1.0)
    end

    @testset "Movement Strategies" begin
        N = 10
        world_size = 100.0f0
        
        for strategy in [EcoEvolutionSim.RANDOM_WALK, EcoEvolutionSim.LANGEVIN, 
                         EcoEvolutionSim.CORRELATED_RW, EcoEvolutionSim.ACTIVE_BROWNIAN]
            
            config = Config(n_agents=N, world_size=world_size, movement_strategy=strategy)
            sim = init_simulation(config)
            
            # Initial positions
            x0 = copy(sim.agents.x)
            y0 = copy(sim.agents.y)
            
            EcoEvolutionSim.movement_step!(sim)
            
            # Check if agents moved
            @test any(sim.agents.x .!= x0) || any(sim.agents.y .!= y0)
            
            # Check boundaries
            @test all(0.0 .<= sim.agents.x .<= world_size)
            @test all(0.0 .<= sim.agents.y .<= world_size)
        end
    end

    @testset "Reflective Boundaries" begin
        # Specifically test reflection logic
        N = 1
        world_size = 10.0f0
        config = Config(n_agents=N, world_size=world_size, movement_strategy=EcoEvolutionSim.LANGEVIN, noise_strength=5.0f0)
        sim = init_simulation(config)
        
        # Force agent to boundary
        sim.agents.x[1] = 0.1f0
        sim.agents.vx[1] = -1.0f0 # Moving left
        
        EcoEvolutionSim.movement_step!(sim)
        
        # Should be reflected
        @test sim.agents.x[1] > 0.1f0
        @test sim.agents.vx[1] > 0.0f0
    end

    @testset "Spatial Grid Logic" begin
        # Test cell_index
        @test cell_index(2.5f0, 2.5f0, 5.0f0, 10) == 1
        @test cell_index(7.5f0, 2.5f0, 5.0f0, 10) == 2
        @test cell_index(2.5f0, 7.5f0, 5.0f0, 10) == 11
        
        nx, ny = 10, 10
        table = EcoEvolutionSim.build_neighbor_table(nx, ny)
        c = 5 + 10*(5-1) # 45
        neighbors = table[c]
        @test length(neighbors) == 9
        @test c in neighbors
    end

    @testset "Interaction Mechanics" begin
        N = 2
        world_size = 10.0f0
        cell_size = 10.0f0
        config = Config(n_agents=N, world_size=world_size, cell_size=cell_size)
        sim = init_simulation(config)
        
        sim.agents.x[1], sim.agents.y[1] = 5.0f0, 5.0f0
        sim.agents.x[2], sim.agents.y[2] = 5.1f0, 5.1f0
        
        # Competition + Mating
        sim.agents.trait[1] = 0.5f0
        sim.agents.trait[2] = 0.5f0
        sim.agents.energy .= 1.0f0
        
        EcoEvolutionSim.build_cell_grid!(sim)
        EcoEvolutionSim.compute_interactions!(sim)
        
        @test sim.agents.energy[1] > 1.0f0 # Net increase
        
        # Predation
        sim.agents.trait[1] = 0.8f0
        sim.agents.trait[2] = 0.2f0
        sim.agents.energy .= 1.0f0
        
        EcoEvolutionSim.build_cell_grid!(sim)
        EcoEvolutionSim.compute_interactions!(sim)
        
        @test sim.agents.energy[1] > 1.0f0
        @test sim.agents.energy[2] < 1.0f0
    end

end
