@testset "Configuration & Core Structures" begin
    @testset "Config Parsing" begin
        # Create a temporary config file
        temp_config = tempname() * ".toml"
        open(temp_config, "w") do io
            write(io, """
            n_agents = 1000
            world_size = 100.0
            strategy = "LANGEVIN"
            """)
        end
        
        cfg = load_config(temp_config)
        @test cfg.n_agents == 1000
        @test cfg.world_size == 100.0f0
        @test cfg.strategy == LANGEVIN
        
        rm(temp_config)
    end
    
    @testset "Enum Conversion" begin
        @test convert(MovementStrategy, "RANDOM_WALK") == RANDOM_WALK
        @test convert(MovementStrategy, "LANGEVIN") == LANGEVIN
        @test_throws ErrorException convert(MovementStrategy, "INVALID")
    end
    
    @testset "Agent Initialization" begin
        traits = (traitA = rand(Float32, 100),)
        agents = Agents(100, 50.0f0, traits)
        
        @test length(agents.x) == 100
        @test length(agents.y) == 100
        @test length(agents.gender) == 100
        @test length(agents.energy) == 100
        @test length(agents.age) == 100
        @test length(agents.alive) == 100
        
        @test all(0.0f0 .<= agents.x .<= 50.0f0)
        @test all(0.0f0 .<= agents.y .<= 50.0f0)
        @test all(in.(agents.gender, Ref((UInt8(0), UInt8(1)))))
        @test all(agents.alive)
        @test agents.traits === traits
    end
end