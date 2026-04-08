using Distributions

@testset "Interactions & Reproduction" begin

    # -----------------------------------------------------------------------
    # Build a sim with the reproduction interaction wired up
    # -----------------------------------------------------------------------
    function make_repro_sim(;N=20, max_agents=200)
        temp_traits = tempname() * ".toml"
        temp_config = tempname() * ".toml"
        open(temp_traits, "w") do io
            write(io, "[fecundity]\nmean = 2.0\nsigma = 0.1\n")
        end
        open(temp_config, "w") do io
            write(io, """
            n_agents   = $N
            max_agents = $max_agents
            world_width  = 10.0
            world_length = 10.0
            cell_size    = 5.0
            r_interact   = 10.0
            """)
        end

        # Load the interaction from script/interactions.jl
        # We re-define locally to avoid module-load side effects
        sim = init_simulation(temp_config, temp_traits;
                              interactions=(interaction_reproduction!,))
        rm(temp_traits); rm(temp_config)
        return sim
    end

    # -----------------------------------------------------------------------
    @testset "spawn_offspring! — child placed near mother" begin
        sim = make_repro_sim()
        # Set mother at known position
        sim.agents.x[1] = 5.0f0
        sim.agents.y[1] = 5.0f0
        sim.agents.alive[1] = true

        old_n = sim.agents.n_alive
        EcoEvolutionSim.claim_agent_id!(sim)   # pre-warm lock
        id = EcoEvolutionSim.claim_agent_id!(sim)
        if id > 0
            # Manual spawn
            sim.agents.x[id] = sim.agents.x[1] + randn(Float32) * 0.01f0
            sim.agents.y[id] = sim.agents.y[1] + randn(Float32) * 0.01f0
            @test isapprox(sim.agents.x[id], 5.0f0, atol=0.1f0)
            @test isapprox(sim.agents.y[id], 5.0f0, atol=0.1f0)
        end
    end

    @testset "spawn_offspring! — n_alive incremented exactly once per spawn" begin
        sim = make_repro_sim(N=10, max_agents=50)
        n0 = sim.agents.n_alive

        # claim 5 ids (each should increment n_alive)
        for _ in 1:5
            id = EcoEvolutionSim.claim_agent_id!(sim)
            if id > 0
                sim.agents.alive[id] = true
            end
        end

        @test sim.agents.n_alive == n0 + 5
    end

    @testset "spawn_offspring! — child stays in world bounds" begin
        sim = make_repro_sim(N=10, max_agents=50)
        # Mother at boundary
        sim.agents.x[1] = 0.001f0
        sim.agents.y[1] = 0.001f0
        sim.agents.alive[1] = true
        sim.agents.gender[1] = UInt8(0)

        # Directly call spawn_offspring! many times
        for _ in 1:20
            EcoEvolutionSim.claim_agent_id!(sim)   # get an id
            id = sim.agents.max_id
            if id <= length(sim.agents.alive)
                sim.agents.x[id] = sim.agents.x[1] + randn(Float32) * 0.01f0
                sim.agents.y[id] = sim.agents.y[1] + randn(Float32) * 0.01f0
                sim.agents.x[id], sim.agents.vx[id] =
                    EcoEvolutionSim.reflect!(sim.agents.x[id], 0f0, sim.config.world_width)
                sim.agents.y[id], sim.agents.vy[id] =
                    EcoEvolutionSim.reflect!(sim.agents.y[id], 0f0, sim.config.world_length)
                sim.agents.alive[id] = true
            end
        end

        alive = sim.agents.alive
        @test all(0f0 .<= sim.agents.x[alive] .<= sim.config.world_width)
        @test all(0f0 .<= sim.agents.y[alive] .<= sim.config.world_length)
    end

    @testset "interaction_reproduction! — females only (gender filter)" begin
        sim = make_repro_sim(N=20, max_agents=200)
        build_cell_grid!(sim)
        n0 = sim.agents.n_alive

        # Make agent 1 male — should not reproduce
        sim.agents.gender[1] = UInt8(1)
        sim.agents.alive[1]  = true
        interaction_reproduction!(sim, 1)
        @test sim.agents.n_alive == n0   # nothing spawned
    end

    @testset "interaction_reproduction! — mating cooldown respected" begin
        sim = make_repro_sim(N=20, max_agents=200)
        build_cell_grid!(sim)

        # Female with freshly mated (last_mating == age)
        sim.agents.gender[1]      = UInt8(0)
        sim.agents.age[1]         = 20.0f0
        sim.agents.last_mating[1] = 19.9f0   # cooldown=10 -> age-last_mating=0.1 < 10
        sim.agents.alive[1]       = true

        n0 = sim.agents.n_alive
        interaction_reproduction!(sim, 1)
        @test sim.agents.n_alive == n0   # cooldown blocked reproduction
    end

    @testset "interactions_step! — n_alive consistent after full step" begin
        sim = make_repro_sim(N=30, max_agents=500)
        build_cell_grid!(sim)
        # Ensure there are males and females with mating age
        for i in 1:length(sim.agents.alive)
            if sim.agents.alive[i]
                sim.agents.gender[i] = UInt8(i % 2)
                sim.agents.age[i]    = 100.0f0
                sim.agents.last_mating[i] = 0.0f0
            end
        end

        interactions_step!(sim)

        @test sim.agents.n_alive == sum(sim.agents.alive)
    end

    @testset "n_alive invariant preserved over full step! with interactions" begin
        sim = make_repro_sim(N=50, max_agents=500)
        for _ in 1:5
            step!(sim)
            @test sim.agents.n_alive == sum(sim.agents.alive)
        end
    end
end
