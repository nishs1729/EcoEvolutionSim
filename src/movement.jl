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
        world = config.world_size
        speed = config.base_speed

        # Pre-calculated constant
        half = 0.5f0 * speed

        x, y = agents.x, agents.y
        vx, vy = agents.vx, agents.vy
        alive = agents.alive

        Threads.@threads for i in eachindex(x)
            @inbounds if alive[i]
                xi, yi = x[i] + rand(Float32) * speed - half, y[i] + rand(Float32) * speed - half

                xi, vxi = reflect!(xi, vx[i], world)
                yi, vyi = reflect!(yi, vy[i], world)

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
        world = config.world_size
        noise = config.noise_strength
        friction = config.friction

        damp = 1f0 - friction
        halfnoise = 0.5f0 * noise

        x, y = agents.x, agents.y
        vx, vy = agents.vx, agents.vy
        alive = agents.alive

        Threads.@threads for i in eachindex(x)
            @inbounds if alive[i]
                vxi = vx[i] * damp + noise * rand(Float32) - halfnoise
                vyi = vy[i] * damp + noise * rand(Float32) - halfnoise

                xi, vxi = reflect!(x[i] + vxi, vxi, world)
                yi, vyi = reflect!(y[i] + vyi, vyi, world)

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
        world = config.world_size
        noise = config.noise_strength
        speed = config.base_speed

        halfnoise = 0.5f0 * noise
        π32 = Float32(π)
        twopi = 2f0 * π32

        x, y = agents.x, agents.y
        theta = agents.theta
        alive = agents.alive

        Threads.@threads for i in eachindex(x)
            @inbounds if alive[i]
                θ = theta[i] + noise * rand(Float32) - halfnoise
                s, c = sincos(θ)
                xi, yi = x[i] + speed * c, y[i] + speed * s

                if xi < 0f0 || xi > world
                    θ = π32 - θ
                end
                if yi < 0f0 || yi > world
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
        world = config.world_size
        noise = config.noise_strength
        speed = config.base_speed

        π32 = Float32(π)
        twopi = 2f0 * π32

        x, y = agents.x, agents.y
        theta = agents.theta
        alive = agents.alive

        Threads.@threads for i in eachindex(x)
            @inbounds if alive[i]
                θ = theta[i] + noise * randn(Float32)
                s, c = sincos(θ)
                xi, yi = x[i] + speed * c, y[i] + speed * s

                if xi < 0f0 || xi > world
                    θ = π32 - θ
                end
                if yi < 0f0 || yi > world
                    θ = -θ
                end

                theta[i] = mod(θ, twopi)
                x[i], y[i] = xi, yi
            end
        end
    end
end