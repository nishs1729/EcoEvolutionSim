function Agents(max_agents, initial_agents, world_size, traits)
    alive = fill(false, max_agents)
    alive[1:initial_agents] .= true

    n_threads = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads()
    death_buffer = [Int[] for _ in 1:n_threads]

    Agents(
        rand(Float32, max_agents) .* world_size,
        rand(Float32, max_agents) .* world_size,
        zeros(Float32, max_agents),
        zeros(Float32, max_agents),
        rand(UInt8(0):UInt8(1), max_agents),
        rand(Float32, max_agents),
        rand(Float32, max_agents) .* 10f0,   # initial age
        alive,
        traits,
        initial_agents,
        Int[],
        death_buffer
    )
end

function ecology_step!(sim::Simulation)
    @timeit to "kernel_ecology" begin
        agents = sim.agents
        age = agents.age
        alive = agents.alive

        μ0 = 0.0001f0
        k = 0.0001f0

        Threads.@threads :static for i in 1:agents.max_id
            @inbounds if alive[i]
                # Age step
                a = age[i] + 0.01f0
                age[i] = a

                # Death step
                @fastmath mortality = μ0 * exp(k * a)
                if rand(Float32) < mortality
                    alive[i] = false
                    push!(agents.death_buffer[Threads.threadid()], i)
                end
            end
        end

        # Merge death buffers
        for t in 1:length(agents.death_buffer)
            if !isempty(agents.death_buffer[t])
                append!(agents.free_indices, agents.death_buffer[t])
                empty!(agents.death_buffer[t])
            end
        end
    end
end