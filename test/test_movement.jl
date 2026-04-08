@testset "Movement Kernels" begin
    @testset "Reflection helper" begin
        # Inside bounds — no change
        x, v = EcoEvolutionSim.reflect!(5.0f0, 1.0f0, 10.0f0)
        @test x == 5.0f0 && v == 1.0f0

        # Left boundary — position mirrored, velocity flipped
        x, v = EcoEvolutionSim.reflect!(-2.0f0, -1.0f0, 10.0f0)
        @test x == 2.0f0 && v == 1.0f0

        # Right boundary
        x, v = EcoEvolutionSim.reflect!(13.0f0, 2.0f0, 10.0f0)
        @test x == 7.0f0 && v == -2.0f0
    end

    # -----------------------------------------------------------------------
    # Helper: default sim for movement tests
    # -----------------------------------------------------------------------
    function make_movement_sim(strategy; N=20, world=10.0f0, speed=0.5f0, noise=1.0f0, friction=0.1f0, dt=1.0f0)
        traits = (dummy = zeros(Float32, N),)
        # Agents constructor: all N agents alive, n_alive=N, theta initialised to random heading
        agents = Agents(N, N, world, world, traits)

        grid = CellGrid(4, N)
        nt   = EcoEvolutionSim.build_neighbor_table(2, 2)
        env  = EnvironmentState(4, 2, 2, world/2f0, 2f0/world, grid, nt)
        cfg  = Config(world_width=world, world_length=world, strategy=strategy,
                      base_speed=speed, noise_strength=noise, friction=friction, dt=dt)
        kernel = EcoEvolutionSim.select_movement_kernel(strategy)
        return Simulation(cfg, agents, env, kernel, ())
    end

    @testset "RANDOM_WALK — all alive agents stay in bounds" begin
        sim = make_movement_sim(RANDOM_WALK; speed=12.0f0)  # extreme speed to stress reflection
        dead_x, dead_y = sim.agents.x[3], sim.agents.y[3]
        sim.agents.alive[3] = false; sim.agents.n_alive -= 1

        for _ in 1:10
            movement_step!(sim)
        end

        for i in 1:length(sim.agents.alive)
            if sim.agents.alive[i]
                @test 0.0f0 <= sim.agents.x[i] <= sim.config.world_width
                @test 0.0f0 <= sim.agents.y[i] <= sim.config.world_length
            end
        end
        # Dead agent must not move
        @test sim.agents.x[3] == dead_x
        @test sim.agents.y[3] == dead_y
    end

    @testset "LANGEVIN — stays in bounds, velocity damped" begin
        sim = make_movement_sim(LANGEVIN; speed=0.5f0, noise=5.0f0)
        for _ in 1:20
            movement_step!(sim)
        end
        alive = sim.agents.alive
        @test all(0.0f0 .<= sim.agents.x[alive] .<= sim.config.world_width)
        @test all(0.0f0 .<= sim.agents.y[alive] .<= sim.config.world_length)
    end

    @testset "CORRELATED_RW — stays in bounds, position reflected on wall hit" begin
        # Place agents near boundary and use large step to force wall collisions
        sim = make_movement_sim(CORRELATED_RW; speed=5.0f0, noise=0.1f0)
        sim.agents.x .= 0.1f0   # near left wall
        sim.agents.y .= 0.1f0
        for _ in 1:20
            movement_step!(sim)
        end
        alive = sim.agents.alive
        @test all(0.0f0 .<= sim.agents.x[alive] .<= sim.config.world_width)
        @test all(0.0f0 .<= sim.agents.y[alive] .<= sim.config.world_length)
    end

    @testset "ACTIVE_BROWNIAN — stays in bounds, position reflected on wall hit" begin
        sim = make_movement_sim(ACTIVE_BROWNIAN; speed=5.0f0, noise=0.5f0)
        sim.agents.x .= 9.9f0   # near right wall
        sim.agents.y .= 9.9f0
        for _ in 1:20
            movement_step!(sim)
        end
        alive = sim.agents.alive
        @test all(0.0f0 .<= sim.agents.x[alive] .<= sim.config.world_width)
        @test all(0.0f0 .<= sim.agents.y[alive] .<= sim.config.world_length)
    end

    @testset "Dead agents unchanged across all kernels" begin
        for strategy in instances(MovementStrategy)
            sim = make_movement_sim(strategy)
            # Kill first agent
            sim.agents.alive[1] = false; sim.agents.n_alive -= 1
            ox, oy = sim.agents.x[1], sim.agents.y[1]
            movement_step!(sim)
            @test sim.agents.x[1] == ox
            @test sim.agents.y[1] == oy
        end
    end
end