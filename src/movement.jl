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
        return 2f0*world - x, -v
    end
    return x, v
end

function movement_random_walk!(sim::Simulation)
    agents = sim.agents
    config = sim.config
    world = config.world_size
    speed = config.base_speed
    half = 0.5f0 * speed

    x = agents.x
    y = agents.y
    vx = agents.vx
    vy = agents.vy
    alive = agents.alive

    if length(x) > MIN_PARALLEL_N
        Threads.@threads for i in eachindex(x)
            if !alive[i]
                continue
            end
            rng = Random.default_rng()
            @inbounds begin
                xi, yi = x[i], y[i]
                vxi, vyi = vx[i], vy[i]

                xi += rand(rng, Float32)*speed - half
                yi += rand(rng, Float32)*speed - half

                xi, vxi = reflect!(xi, vxi, world)
                yi, vyi = reflect!(yi, vyi, world)

                x[i], y[i] = xi, yi
                vx[i], vy[i] = vxi, vyi
            end
        end
    else
        rng = Random.default_rng()
        for i in eachindex(x)
            if !alive[i]
                continue
            end
            @inbounds begin
                xi, yi = x[i], y[i]
                vxi, vyi = vx[i], vy[i]

                xi += rand(rng, Float32)*speed - half
                yi += rand(rng, Float32)*speed - half

                xi, vxi = reflect!(xi, vxi, world)
                yi, vyi = reflect!(yi, vyi, world)

                x[i], y[i] = xi, yi
                vx[i], vy[i] = vxi, vyi
            end
        end
    end
end

function movement_langevin!(sim::Simulation)
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

    if length(x) > MIN_PARALLEL_N
        Threads.@threads for i in eachindex(x)
            if !alive[i]
                continue
            end
            rng = Random.default_rng()
            @inbounds begin
                xi, yi = x[i], y[i]
                vxi, vyi = vx[i], vy[i]

                vxi = vxi*damp + noise*rand(rng, Float32) - halfnoise
                vyi = vyi*damp + noise*rand(rng, Float32) - halfnoise

                xi += vxi
                yi += vyi

                xi, vxi = reflect!(xi, vxi, world)
                yi, vyi = reflect!(yi, vyi, world)

                x[i], y[i] = xi, yi
                vx[i], vy[i] = vxi, vyi
            end
        end
    else
        rng = Random.default_rng()
        for i in eachindex(x)
            if !alive[i]
                continue
            end
            @inbounds begin
                xi, yi = x[i], y[i]
                vxi, vyi = vx[i], vy[i]

                vxi = vxi*damp + noise*rand(rng, Float32) - halfnoise
                vyi = vyi*damp + noise*rand(rng, Float32) - halfnoise

                xi += vxi
                yi += vyi

                xi, vxi = reflect!(xi, vxi, world)
                yi, vyi = reflect!(yi, vyi, world)

                x[i], y[i] = xi, yi
                vx[i], vy[i] = vxi, vyi
            end
        end
    end
end

function movement_correlated_rw!(sim::Simulation)
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

    if length(x) > MIN_PARALLEL_N
        Threads.@threads for i in eachindex(x)
            if !alive[i]
                continue
            end
            rng = Random.default_rng()
            @inbounds begin
                xi, yi = x[i], y[i]
                θ = theta[i]

                θ += noise*rand(rng, Float32) - halfnoise
                s, c = sincos(θ)
                xi += speed*c
                yi += speed*s

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
    else
        rng = Random.default_rng()
        for i in eachindex(x)
            if !alive[i]
                continue
            end
            @inbounds begin
                xi, yi = x[i], y[i]
                θ = theta[i]

                θ += noise*rand(rng, Float32) - halfnoise
                s, c = sincos(θ)
                xi += speed*c
                yi += speed*s

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

    if length(x) > MIN_PARALLEL_N
        Threads.@threads for i in eachindex(x)
            if !alive[i]
                continue
            end
            rng = Random.default_rng()
            @inbounds begin
                xi, yi = x[i], y[i]
                θ = theta[i]

                θ += noise * randn(rng, Float32)
                s, c = sincos(θ)
                xi += speed*c
                yi += speed*s

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
    else
        rng = Random.default_rng()
        for i in eachindex(x)
            if !alive[i]
                continue
            end
            @inbounds begin
                xi, yi = x[i], y[i]
                θ = theta[i]

                θ += noise * randn(rng, Float32)
                s, c = sincos(θ)
                xi += speed*c
                yi += speed*s

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