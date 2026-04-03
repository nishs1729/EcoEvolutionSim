@testset "Ecology & Simulation Loop" begin
    @testset "Ecology Step" begin
        N = 5
        traits = (dummy = zeros(Float32, N),)
        agents = Agents(N, 10.0f0, traits)
        
        agents.age .= [1.0f0, 2.0f0, 3.0f0, 4.0f0, 5.0f0]
        agents.alive .= [true, true, false, true, true]
        
        config = Config()
        neighbor_table = NeighborTable(Matrix{Int32}(undef, 0, 0), Int8[])
        env = EnvironmentState(1, 1, 1, 10.0f0, 0.1f0, CellGrid(1, N), neighbor_table)
        sim = Simulation(config, agents, env, x -> nothing)
        
        ecology_step!(sim)
        
        # Check age increments for alive agents
        @test isapprox(agents.age[1], 1.01f0)
        @test isapprox(agents.age[2], 2.01f0)
        @test agents.age[3] == 3.0f0 # dead agent age should not increase
        
        # Note: current ecology_step! has hardcoded μ0 = 0.0, so mortality is 0
        @test all(agents.alive .== [true, true, false, true, true])
    end
    
    @testset "Initialization & Loop" begin
        # Mock traits file
        temp_traits = tempname() * ".toml"
        open(temp_traits, "w") do io
            write(io, """
            [resource_pref]
            mean = 0.5
            sigma = 0.1
            """)
        end
        
        # We can just use default Config via `init_simulation`
        # but it defaults to config.toml, so let's mock one
        temp_config = tempname() * ".toml"
        open(temp_config, "w") do io
            write(io, """
            n_agents = 100
            world_size = 50.0
            cell_size = 5.0
            """)
        end
        
        sim = init_simulation(temp_config, temp_traits)
        
        @test length(sim.agents.x) == 100
        @test sim.env.nx == 10
        @test sim.env.ny == 10
        @test sim.env.ncells == 100
        
        # Run 5 steps
        for _ in 1:5
            step!(sim)
        end
        
        # Verify it didn't crash and grid state seems consistent
        @test sum(sim.env.grid.cell_count) <= 100 # alive agents only
        
        rm(temp_traits)
        rm(temp_config)
    end
end