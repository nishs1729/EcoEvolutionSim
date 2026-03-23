module Fitness

using ..EcoEvolutionSim: Simulation, Config
using ..FitnessRegistry: register_fitness

export compute_interactions!

# Utility function used by multiple modules
@inline function cell_index(x, y, cell_size, nx)
    ix = Int(floor(x / cell_size)) + 1
    iy = Int(floor(y / cell_size)) + 1
    return ix + nx * (iy - 1)
end

@inline function default_interaction_kernel!(i, j, sim::Simulation)
    agents = sim.agents
    config = sim.config

    ti = agents.trait[i]
    tj = agents.trait[j]

    d = ti - tj
    ad = abs(d)

    # Competition
    comp = exp(-(ad^2) / (2f0 * config.competition_sigma^2))
    agents.energy[i] -= 0.01f0 * comp
    agents.energy[j] -= 0.01f0 * comp

    # Mating
    mate = exp(-(ad^2) / (2f0 * config.mating_sigma^2))
    if mate > 0.8f0
        agents.energy[i] += 0.05f0
        agents.energy[j] += 0.05f0
    end

    # Predation
    if d > config.predation_threshold
        agents.energy[i] += 0.1f0
        agents.energy[j] -= 0.1f0
    elseif d < -config.predation_threshold
        agents.energy[j] += 0.1f0
        agents.energy[i] -= 0.1f0
    end
end

# Phase 6: Register the existing fitness as a plugin
# For now, the constructor just returns the kernel function
function default_interaction_constructor(config::Config)
    return default_interaction_kernel!
end

# This needs to be called somewhere, or we can use a __init__ function
function __init__()
    register_fitness("default_interaction", default_interaction_constructor)
end

function compute_interactions!(sim::Simulation, kernel! = default_interaction_kernel!)
    agents = sim.agents
    env = sim.env
    config = sim.config
    N = length(agents.x)

    Threads.@threads for i in 1:N
        xi = agents.x[i]
        yi = agents.y[i]
        c = cell_index(xi, yi, env.cell_size, env.nx)

        for nc in env.neighbor_cells[c]
            start = env.grid.cell_start[nc]
            stop = start + env.grid.cell_count[nc][] - 1
            for k in start:stop
                j = env.grid.agent_index[k]
                if j ≤ i continue end

                dx = agents.x[j] - xi
                dy = agents.y[j] - yi
                r2 = dx * dx + dy * dy

                if r2 < config.interaction_radius2
                    kernel!(i, j, sim)
                end
            end
        end
    end
end

end
