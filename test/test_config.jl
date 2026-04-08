@testset "Configuration & Core Structures" begin
    @testset "Config Parsing" begin
        temp_config = tempname() * ".toml"
        open(temp_config, "w") do io
            write(io, """
            n_agents = 1000
            world_width = 100.0
            world_length = 100.0
            strategy = "LANGEVIN"
            """)
        end

        cfg = load_config(temp_config)
        @test cfg.n_agents == 1000
        @test cfg.world_width == 100.0f0
        @test cfg.world_length == 100.0f0
        @test cfg.strategy == LANGEVIN

        rm(temp_config)
    end

    @testset "Enum Conversion — all strategies" begin
        @test convert(MovementStrategy, "RANDOM_WALK")     == RANDOM_WALK
        @test convert(MovementStrategy, "LANGEVIN")        == LANGEVIN
        @test convert(MovementStrategy, "CORRELATED_RW")   == CORRELATED_RW
        @test convert(MovementStrategy, "ACTIVE_BROWNIAN") == ACTIVE_BROWNIAN
        # Mixed/lower case is accepted (uppercased internally)
        @test convert(MovementStrategy, "random_walk")     == RANDOM_WALK
        @test convert(MovementStrategy, "Langevin")        == LANGEVIN
        # Truly invalid string must throw
        @test_throws ErrorException convert(MovementStrategy, "INVALID")
        @test_throws ErrorException convert(MovementStrategy, "WALK")
    end

    @testset "Agent Initialization — shapes & n_alive" begin
        traits = (traitA = rand(Float32, 100),)
        agents = Agents(100, 60, 50.0f0, 50.0f0, traits)

        @test length(agents.x)      == 100
        @test length(agents.y)      == 100
        @test length(agents.gender) == 100
        @test length(agents.energy) == 100
        @test length(agents.age)    == 100
        @test length(agents.alive)  == 100

        @test all(0.0f0 .<= agents.x .<= 50.0f0)
        @test all(0.0f0 .<= agents.y .<= 50.0f0)
        @test all(in.(agents.gender, Ref((UInt8(0), UInt8(1)))))

        # First 60 alive; rest dead
        @test sum(agents.alive) == 60
        @test agents.max_id     == 60
        # n_alive counter must match
        @test agents.n_alive == 60
        @test agents.n_alive == sum(agents.alive)

        @test agents.traits === traits
    end
end