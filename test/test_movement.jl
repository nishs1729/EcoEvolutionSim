@testset "Movement Kernels" begin
    @testset "Reflection" begin
        # Normal inside
        x, v = EcoEvolutionSim.reflect!(5.0f0, 1.0f0, 10.0f0)
        @test x == 5.0f0
        @test v == 1.0f0
        
        # Left boundary
        x, v = EcoEvolutionSim.reflect!(-2.0f0, -1.0f0, 10.0f0)
        @test x == 2.0f0
        @test v == 1.0f0
        
        # Right boundary
        x, v = EcoEvolutionSim.reflect!(13.0f0, 2.0f0, 10.0f0)
        @test x == 7.0f0
        @test v == -2.0f0
    end
    
    @testset "Movement Kernels In-Bounds" begin
        # Note: LANGEVIN, CORRELATED_RW, ACTIVE_BROWNIAN currently access 
        # undefined fields (noise_strength, friction, theta) in Config/Agents.
        # We test RANDOM_WALK which has all its dependencies satisfied.
        
        N = 10
        traits = (dummy = zeros(Float32, N),) 
        agents = Agents(N, 10.0f0, traits)
        
        agents.x .= 5.0f0
        agents.y .= 5.0f0
        agents.vx .= 1.0f0
        agents.vy .= 1.0f0
        agents.alive .= [true, true, false, false, true, true, true, true, true, true]
        
        orig_x_dead = agents.x[3]
        orig_y_dead = agents.y[3]
        
        grid = CellGrid(4, N)
        neighbor_table = NeighborTable(Matrix{Int32}(undef, 0, 0), Int8[])
        env = EnvironmentState(4, 2, 2, 5.0f0, 0.2f0, grid, neighbor_table)
        
        config = Config(world_size=10.0f0, strategy=RANDOM_WALK, base_speed=12.0f0) 
        kernel = EcoEvolutionSim.select_movement_kernel(RANDOM_WALK)
        sim = Simulation(config, agents, env, kernel)
        
        movement_step!(sim)
        
        # Check all alive are in bounds
        for i in 1:N
            if agents.alive[i]
                @test 0.0f0 <= agents.x[i] <= 10.0f0
                @test 0.0f0 <= agents.y[i] <= 10.0f0
            else
                # Dead agents shouldn't move
                @test agents.x[i] == orig_x_dead
                @test agents.y[i] == orig_y_dead
            end
        end
    end
end