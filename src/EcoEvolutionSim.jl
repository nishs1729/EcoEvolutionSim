module EcoEvolutionSim

include("simulator.jl")

export init_simulation, step!, load_config
export TRAIT_SPECS
export Simulation, Agents, CellGrid, EnvironmentState
export cell_index
export Config, MovementStrategy, RANDOM_WALK, LANGEVIN, CORRELATED_RW, ACTIVE_BROWNIAN

end
