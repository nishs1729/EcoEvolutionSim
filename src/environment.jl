module Environment

using ..EcoEvolutionSim: Config, CellGrid

export EnvironmentState

struct EnvironmentState
    # Parameters from config that define the environment
    world_size::Float32
    cell_size::Float32
    
    # Grid and neighbor table are technically part of the spatial environment
    grid::CellGrid
    neighbor_cells::Vector{Vector{Int}}
    
    # Grid dimensions
    nx::Int
    ny::Int
    ncells::Int
end

function EnvironmentState(config::Config, grid::CellGrid, neighbor_cells::Vector{Vector{Int}}, nx::Int, ny::Int, ncells::Int)
    return EnvironmentState(
        config.sim.world_size,
        config.sim.cell_size,
        grid,
        neighbor_cells,
        nx,
        ny,
        ncells
    )
end

end
