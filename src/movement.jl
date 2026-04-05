function select_movement_kernel(strategy::MovementStrategy)
    if strategy == RANDOM_WALK
        return movement_random_walk!
    elseif strategy == LANGEVIN
        return movement_langevin!
    elseif strategy == CORRELATED_RW
        return movement_correlated_rw!
    elseif strategy == ACTIVE_BROWNIAN
        return movement_active_brownian!
    else
        error("Unknown movement strategy: $strategy")
    end
end

function movement_step!(sim::Simulation)
    sim.movement_kernel!(sim)
end

@inline function reflect!(x, v, world)
    if x < 0f0
        return -x, -v
    elseif x > world
        return 2f0 * world - x, -v
    end
    return x, v
end

function movement_random_walk!(sim::Simulation)
    @timeit to "kernel_rw" begin
        agents = sim.agents
        config = sim.config
        world_width = config.world_width
        world_length = config.world_length
        speed = config.base_speed
        dt = config.dt

        # Pre-calculated constant
        step_dist = speed * dt
        half = 0.5f0 * step_dist

        x, y = agents.x, agents.y
        vx, vy = agents.vx, agents.vy
        alive = agents.alive

        Threads.@threads for i in 1:agents.max_id
            @inbounds if alive[i]
                xi, yi = x[i] + rand(Float32) * step_dist - half, y[i] + rand(Float32) * step_dist - half

                xi, vxi = reflect!(xi, vx[i], world_width)
                yi, vyi = reflect!(yi, vy[i], world_length)

                x[i], y[i] = xi, yi
                vx[i], vy[i] = vxi, vyi
            end
        end
    end
end

function movement_langevin!(sim::Simulation)
    @timeit to "kernel_langevin" begin
        agents = sim.agents
        config = sim.config
        world_width = config.world_width
        world_length = config.world_length
        noise = config.noise_strength
        friction = config.friction
        dt = config.dt

        damp = 1f0 - friction * dt
        step_noise = noise * sqrt(dt)
        halfnoise = 0.5f0 * step_noise

        x, y = agents.x, agents.y
        vx, vy = agents.vx, agents.vy
        alive = agents.alive

        Threads.@threads for i in 1:agents.max_id
            @inbounds if alive[i]
                vxi = vx[i] * damp + step_noise * rand(Float32) - halfnoise
                vyi = vy[i] * damp + step_noise * rand(Float32) - halfnoise

                xi, vxi = reflect!(x[i] + vxi * dt, vxi, world_width)
                yi, vyi = reflect!(y[i] + vyi * dt, vyi, world_length)

                x[i], y[i] = xi, yi
                vx[i], vy[i] = vxi, vyi
            end
        end
    end
end

function movement_correlated_rw!(sim::Simulation)
    @timeit to "kernel_crw" begin
        agents = sim.agents
        config = sim.config
        world_width = config.world_width
        world_length = config.world_length
        noise = config.noise_strength
        speed = config.base_speed
        dt = config.dt

        step_dist = speed * dt
        step_noise = noise * sqrt(dt)
        halfnoise = 0.5f0 * step_noise
        π32 = Float32(π)
        twopi = 2f0 * π32

        x, y = agents.x, agents.y
        theta = agents.theta
        alive = agents.alive

        Threads.@threads for i in 1:agents.max_id
            @inbounds if alive[i]
                θ = theta[i] + step_noise * rand(Float32) - halfnoise
                s, c = sincos(θ)
                xi, yi = x[i] + step_dist * c, y[i] + step_dist * s

                if xi < 0f0
                    xi = -xi
                    θ = π32 - θ
                elseif xi > world_width
                    xi = 2f0 * world_width - xi
                    θ = π32 - θ
                end
                if yi < 0f0
                    yi = -yi
                    θ = -θ
                elseif yi > world_length
                    yi = 2f0 * world_length - yi
                    θ = -θ
                end

                theta[i] = mod(θ, twopi)
                x[i], y[i] = xi, yi
            end
        end
    end
end

function movement_active_brownian!(sim::Simulation)
    @timeit to "kernel_abm" begin
        agents = sim.agents
        config = sim.config
        world_width = config.world_width
        world_length = config.world_length
        noise = config.noise_strength
        speed = config.base_speed
        dt = config.dt
        
        step_dist = speed * dt
        step_noise = noise * sqrt(dt)

        π32 = Float32(π)
        twopi = 2f0 * π32

        x, y = agents.x, agents.y
        theta = agents.theta
        alive = agents.alive

        Threads.@threads for i in 1:agents.max_id
            @inbounds if alive[i]
                θ = theta[i] + step_noise * randn(Float32)
                s, c = sincos(θ)
                xi, yi = x[i] + step_dist * c, y[i] + step_dist * s

                if xi < 0f0
                    xi = -xi
                    θ = π32 - θ
                elseif xi > world_width
                    xi = 2f0 * world_width - xi
                    θ = π32 - θ
                end
                if yi < 0f0
                    yi = -yi
                    θ = -θ
                elseif yi > world_length
                    yi = 2f0 * world_length - yi
                    θ = -θ
                end

                theta[i] = mod(θ, twopi)
                x[i], y[i] = xi, yi
            end
        end
    end
end