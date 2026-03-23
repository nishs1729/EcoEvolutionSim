module EcoEvolutionSim

include("registry/fitness_registry.jl")
include("registry/trait_registry.jl")
include("simulator.jl")

export init_simulation, step!, load_config, compute_interactions!
export Simulation, Agents, CellGrid, EnvironmentState
export cell_index
export Config, MovementStrategy, RANDOM_WALK, LANGEVIN, CORRELATED_RW, ACTIVE_BROWNIAN
export FitnessRegistry, TraitRegistry

end
