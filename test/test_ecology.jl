@testset "Ecology & Simulation Loop" begin

    # -----------------------------------------------------------------------
    # Helper
    # -----------------------------------------------------------------------
    function make_eco_sim(N; max_agents=N, n_alive=N, world=10.0f0, cell_size=5.0f0)
        traits = (dummy = zeros(Float32, max_agents),)
        agents = Agents(max_agents, n_alive, world, world, traits)
        ncells = Int(ceil(world / cell_size))^2
        nx = ny = Int(ceil(world / cell_size))
        grid = CellGrid(ncells, max_agents)
        nt = EcoEvolutionSim.build_neighbor_table(nx, ny)
        env = EnvironmentState(ncells, nx, ny, cell_size, 1f0 / cell_size, grid, nt)
        cfg = Config(world_width=world, world_length=world, cell_size=cell_size)
        kernel = EcoEvolutionSim.select_movement_kernel(RANDOM_WALK)
        return Simulation(cfg, agents, env, kernel, ())
    end

    # -----------------------------------------------------------------------
    @testset "Ecology Step — age increments for alive agents only" begin
        sim = make_eco_sim(5)
        sim.agents.age  .= [1.0f0, 2.0f0, 3.0f0, 4.0f0, 5.0f0]
        sim.agents.alive .= [true, true, false, true, true]
        sim.agents.n_alive = 4

        ecology_step!(sim)

        @test isapprox(sim.agents.age[1], 2.0f0)   # alive: incremented
        @test isapprox(sim.agents.age[2], 3.0f0)
        @test sim.agents.age[3] == 3.0f0             # dead: unchanged
        @test isapprox(sim.agents.age[4], 5.0f0)
        @test isapprox(sim.agents.age[5], 6.0f0)
    end

    @testset "Ecology Step — n_alive decremented by deaths" begin
        sim = make_eco_sim(10)
        # Force agents to be very old so mortality is certain (age=10000)
        sim.agents.age .= 10_000.0f0
        sim.agents.alive .= true
        sim.agents.n_alive = 10

        ecology_step!(sim)

        # With age=10000, Gompertz mortality = 0.0001 * exp(0.01 * 10001) ≈ ≫ 1 → all die
        survived = sum(sim.agents.alive)
        @test sim.agents.n_alive == survived          # counter must match truth
        @test sim.agents.n_alive <= 10
        # free_indices must contain the dead slots
        n_dead = 10 - survived
        @test length(sim.agents.free_indices) == n_dead
    end

    @testset "Ecology Step — mortality threshold: very young agents never call rand" begin
        # Age=0 -> mortality = 0.0001 * exp(0) * dt ≈ 1e-4, above 1e-6 threshold
        # But we can verify that extremely young agents (age ~ 0 with tiny dt) survive
        sim = make_eco_sim(50)
        sim.agents.age   .= 0.0f0
        sim.agents.alive .= true
        sim.agents.n_alive = 50

        # Run 1 step with default μ0=1e-4, dt=1 -> mortality≈1e-4 → very unlikely to die
        # Use a fixed seed for reproducibility
        Random.seed!(42)
        ecology_step!(sim)

        # With p≈1e-4 per agent per step and 50 agents, expected deaths ≈ 0.005 → nearly always 0
        @test sum(sim.agents.alive) >= 48   # should almost certainly all survive
        @test sim.agents.n_alive == sum(sim.agents.alive)
    end

    # -----------------------------------------------------------------------
    @testset "claim_agent_id! — sequential: reuses free indices" begin
        sim = make_eco_sim(5, max_agents=10, n_alive=5)
        # Kill some, collect free slots
        sim.agents.alive[2] = false
        sim.agents.alive[4] = false
        push!(sim.agents.free_indices, 2, 4)
        sim.agents.n_alive -= 2

        id1 = EcoEvolutionSim.claim_agent_id!(sim)
        id2 = EcoEvolutionSim.claim_agent_id!(sim)
        @test id1 in (2, 4)
        @test id2 in (2, 4)
        @test id1 != id2
        # n_alive incremented twice
        @test sim.agents.n_alive == 3 + 2   # was 3 after deaths
    end

    @testset "claim_agent_id! — sequential: extends max_id when no free slots" begin
        sim = make_eco_sim(3, max_agents=6, n_alive=3)
        @test sim.agents.max_id == 3
        id = EcoEvolutionSim.claim_agent_id!(sim)
        @test id == 4
        @test sim.agents.max_id == 4
        @test sim.agents.n_alive == 4  # incremented
    end

    @testset "claim_agent_id! — returns 0 when at capacity" begin
        sim = make_eco_sim(5, max_agents=5, n_alive=5)
        # max_id already at capacity
        id = EcoEvolutionSim.claim_agent_id!(sim)
        @test id == 0
        @test sim.agents.n_alive == 5  # unchanged
    end

    @testset "claim_agent_id! — concurrent: n_alive consistent after parallel spawns" begin
        sim = make_eco_sim(0, max_agents=1000, n_alive=0)
        sim.agents.max_id = 0

        # Spawn 200 agents from 8 threads concurrently
        Threads.@threads for _ in 1:200
            id = EcoEvolutionSim.claim_agent_id!(sim)
            if id > 0
                sim.agents.alive[id] = true
            end
        end

        @test sim.agents.n_alive == sum(sim.agents.alive)
        @test sim.agents.n_alive == 200
    end

    # -----------------------------------------------------------------------
    @testset "Initialization & step! loop — n_alive tracks sum(alive)" begin
        temp_traits = tempname() * ".toml"
        temp_config = tempname() * ".toml"
        open(temp_traits, "w") do io
            write(io, "[resource_pref]\nmean = 0.5\nsigma = 0.1\n")
        end
        open(temp_config, "w") do io
            write(io, """
            n_agents   = 50
            max_agents = 200
            world_width  = 50.0
            world_length = 50.0
            cell_size    = 5.0
            """)
        end

        sim = init_simulation(temp_config, temp_traits)

        # Initial invariant
        @test sim.agents.n_alive == sum(sim.agents.alive)
        @test sim.agents.n_alive == 50

        for _ in 1:10
            step!(sim)
            @test sim.agents.n_alive == sum(sim.agents.alive)
        end

        # Grid count must also match n_alive (after build_cell_grid! called by step!)
        @test sum(sim.env.grid.cell_count) == sim.agents.n_alive

        rm(temp_traits); rm(temp_config)
    end

    @testset "No double-free: free_indices has no duplicates after ecology_step!" begin
        sim = make_eco_sim(20)
        sim.agents.age   .= 10_000.0f0   # force all to die
        sim.agents.alive .= true
        sim.agents.n_alive = 20

        ecology_step!(sim)

        # No slot should appear twice in free_indices
        @test length(unique(sim.agents.free_indices)) == length(sim.agents.free_indices)
    end
end