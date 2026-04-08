function Agents(max_agents, initial_agents, world_width, world_length, traits)
    alive = fill(false, max_agents)
    alive[1:initial_agents] .= true

    n_threads = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads()
    death_buffer = [Int[] for _ in 1:n_threads]

    Agents(
        rand(Float32, max_agents) .* world_width,
        rand(Float32, max_agents) .* world_length,
        zeros(Float32, max_agents),
        zeros(Float32, max_agents),
        rand(Float32, max_agents) .* (2f0 * Float32(π)),  # theta: random initial heading
        rand(UInt8(0):UInt8(1), max_agents),
        rand(Float32, max_agents),
        rand(Float32, max_agents) .* 10f0,   # initial age
        zeros(Float32, max_agents),          # last_mating
        alive,
        traits,
        initial_agents,
        initial_agents,                      # n_alive starts equal to initial_agents
        Int[],
        death_buffer
    )
end

function ecology_step!(sim::Simulation)
    @timeit to "kernel_ecology" begin
        agents = sim.agents
        age = agents.age
        alive = agents.alive
        energy = agents.energy
        dt = sim.config.dt
        metabolic_drain = sim.config.metabolic_rate * dt

        μ0 = 0.0001f0
        k = 0.01f0
        mortality_threshold = 1f-6  # skip rand() when mortality is negligibly small

        Threads.@threads :static for i in 1:agents.max_id
            @inbounds if alive[i]
                # Age step
                a = age[i] + dt
                age[i] = a

                # Metabolism step
                energy[i] -= metabolic_drain

                # Death step — early-out avoids rand() call for young agents
                @fastmath mortality = μ0 * exp(k * a) * dt
                if mortality > mortality_threshold && rand(Float32) < mortality
                    alive[i] = false
                    push!(agents.death_buffer[Threads.threadid()], i)
                end
            end
        end

        # Merge death buffers and update n_alive counter
        for t in 1:length(agents.death_buffer)
            if !isempty(agents.death_buffer[t])
                agents.n_alive -= length(agents.death_buffer[t])
                append!(agents.free_indices, agents.death_buffer[t])
                empty!(agents.death_buffer[t])
            end
        end
    end
end

const SPAWN_LOCK = Threads.SpinLock()

function claim_agent_id!(sim::Simulation)
    lock(SPAWN_LOCK)
    try
        if !isempty(sim.agents.free_indices)
            id = pop!(sim.agents.free_indices)
            sim.agents.n_alive += 1   # inside lock — thread-safe
            return id
        else
            id = sim.agents.max_id + 1
            if id <= length(sim.agents.alive)
                sim.agents.max_id = id
                sim.agents.n_alive += 1   # inside lock — thread-safe
                return id
            else
                return 0 # Out of bounds
            end
        end
    finally
        unlock(SPAWN_LOCK)
    end
end

function interactions_step!(sim::Simulation)
    @timeit to "kernel_interactions" begin
        if !isempty(sim.interactions)
            Threads.@threads :static for i in 1:sim.agents.max_id
                @inbounds if sim.agents.alive[i]
                    for interact! in sim.interactions
                        interact!(sim, i)
                    end
                end
            end
        end
    end
end