function Agents(N, world_size, traits_dict::Dict{String, Vector{Float32}})
    Agents(
        rand(Float32, N) .* world_size,
        rand(Float32, N) .* world_size,
        zeros(Float32, N),
        zeros(Float32, N),
        rand(UInt8(0):UInt8(1), N),
        rand(Float32, N),
        rand(Float32, N) .* 10f0,   # initial age
        trues(N),                   # alive
        traits_dict
    )
end

function age_step!(sim::Simulation)
    agents = sim.agents
    age = agents.age
    alive = agents.alive

    if length(age) > MIN_PARALLEL_N
        Threads.@threads for i in eachindex(age)
            if alive[i]
                age[i] += 0.01f0
            end
        end
    else
        for i in eachindex(age)
            if alive[i]
                age[i] += 0.01f0
            end
        end
    end
end

function death_step!(sim::Simulation)
    agents = sim.agents
    age = agents.age
    alive = agents.alive

    μ0 = 0.000f0
    k = 0.000f0

    if length(age) > MIN_PARALLEL_N
        Threads.@threads for i in eachindex(age)
            if !alive[i]
                continue
            end

            a = age[i]
            @fastmath mortality = μ0 * exp(k*a)
            if rand(Float32) < mortality
                alive[i] = false
            end
        end
    else
        rng = Random.default_rng()
        for i in eachindex(age)
            if !alive[i]
                continue
            end

            a = age[i]
            @fastmath mortality = μ0 * exp(k*a)
            if rand(rng, Float32) < mortality
                alive[i] = false
            end
        end
    end
end