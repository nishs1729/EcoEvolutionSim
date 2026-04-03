using Distributions
using Random
using EcoEvolutionSim

# Parameters
const p_mate = 0.1f0
const mating_cooldown = 10.0f0 # 2 days

function interaction_reproduction!(sim::Simulation, i::Int)
    agents = sim.agents

    # Check if agent is a female (0)
    if agents.gender[i] != 0
        return
    end

    # Check cooldown
    if agents.age[i] - agents.last_mating[i] < mating_cooldown
        return
    end

    # Search for a mate within r_interact
    r_interact = sim.config.r_interact
    males = Int[]

    for j in nearby_agents(sim, i, r_interact)
        if agents.gender[j] == 1
            push!(males, j)
        end
    end

    if !isempty(males) && rand(Float32) < p_mate
        father = rand(males)
        fecundity = agents.traits.fecundity[i]

        # Determine number of offspring
        n_offspring = rand(Poisson(fecundity))

        # Reset cooldown tracking
        agents.last_mating[i] = agents.age[i]

        # Spawn each offspring sequentially
        for k in 1:n_offspring
            spawn_offspring!(sim, i, father)
        end
    end
end

function spawn_offspring!(sim::Simulation, mother_id::Int, father_id::Int)
    id = EcoEvolutionSim.claim_agent_id!(sim)
    if id == 0
        return # Max capacity reached
    end

    agents = sim.agents

    # Placed near mother with a slight micro-jitter
    agents.x[id] = agents.x[mother_id] + randn(Float32) * 0.01f0
    agents.y[id] = agents.y[mother_id] + randn(Float32) * 0.01f0

    # Bounds check via reflect to prevent escaping
    world_width = sim.config.world_width
    world_length = sim.config.world_length
    agents.x[id], agents.vx[id] = EcoEvolutionSim.reflect!(agents.x[id], 0.0f0, world_width)
    agents.y[id], agents.vy[id] = EcoEvolutionSim.reflect!(agents.y[id], 0.0f0, world_length)

    agents.gender[id] = rand(UInt8(0):UInt8(1))
    agents.energy[id] = 10.0f0
    agents.age[id] = 0.0f0
    agents.last_mating[id] = 0.0f0
    agents.alive[id] = true

    # Inheritance with 10% mutation magnitude randomly drawn from a parent
    for key in keys(agents.traits)
        prop = getproperty(agents.traits, key)
        mother_trait = prop[mother_id]
        father_trait = prop[father_id]
        baseline = rand(Bool) ? mother_trait : father_trait
        prop[id] = max(0.0f0, baseline + randn(Float32) * 0.1f0 * baseline)
    end
end
