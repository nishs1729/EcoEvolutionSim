function Agents(N, world_size, traits)
    Agents(
        rand(Float32, N) .* world_size,
        rand(Float32, N) .* world_size,
        zeros(Float32, N),
        zeros(Float32, N),
        rand(UInt8(0):UInt8(1), N),
        rand(Float32, N),
        rand(Float32, N) .* 10f0,   # initial age
        fill(true, N),                   # alive
        traits
    )
end

function ecology_step!(sim::Simulation)
    @timeit to "kernel_ecology" begin
        agents = sim.agents
        age = agents.age
        alive = agents.alive

        μ0 = 0.000f0
        k = 0.000f0

        # Pre-check if death is even possible
        death_enabled = μ0 > 0f0

        if length(age) > MIN_PARALLEL_N
            Threads.@threads for i in eachindex(age)
                @inbounds if alive[i]
                    # Age step
                    a = age[i] + 0.01f0
                    age[i] = a

                    # Death step
                    if death_enabled
                        @fastmath mortality = μ0 * exp(k*a)
                        if rand(Float32) < mortality
                            alive[i] = false
                        end
                    end
                end
            end
        else
            for i in eachindex(age)
                @inbounds if alive[i]
                    a = age[i] + 0.01f0
                    age[i] = a

                    if death_enabled
                        @fastmath mortality = μ0 * exp(k*a)
                        if rand(Float32) < mortality
                            alive[i] = false
                        end
                    end
                end
            end
        end
    end
end