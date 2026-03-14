############################################################
# Phenotype-space ecological simulation
# Designed for millions of agents
############################################################

using Base.Threads
using Random

############################################################
# Agent storage (Structure of Arrays)
############################################################

struct Agents
    x::Vector{Float32}
    y::Vector{Float32}
    trait::Vector{Float32}
    energy::Vector{Float32}
end

############################################################
# Spatial grid
############################################################

struct CellGrid
    cell_start::Vector{Int}
    cell_count::Vector{Int}
    cell_offset::Vector{Int}
    agent_index::Vector{Int}
end

############################################################
# Simulation object
############################################################

struct Simulation

    agents::Agents

    world_size::Float32
    cell_size::Float32

    nx::Int
    ny::Int
    ncells::Int

    neighbor_cells::Vector{Vector{Int}}

    grid::CellGrid

    interaction_radius2::Float32
    competition_sigma::Float32
    mating_sigma::Float32
    predation_threshold::Float32

end

############################################################
# Convert position → grid cell
############################################################

@inline function cell_index(x,y,cell_size,nx)

    ix = Int(floor(x/cell_size)) + 1
    iy = Int(floor(y/cell_size)) + 1

    return ix + nx*(iy-1)

end

############################################################
# Precompute neighbor cell table
############################################################

function build_neighbor_table(nx,ny)

    ncells = nx*ny
    table = Vector{Vector{Int}}(undef,ncells)

    for c in 1:ncells

        iy = (c-1) ÷ nx + 1
        ix = c - nx*(iy-1)

        list = Int[]

        for dy in -1:1
            for dx in -1:1

                x2 = ix + dx
                y2 = iy + dy

                if 1 ≤ x2 ≤ nx && 1 ≤ y2 ≤ ny
                    push!(list, x2 + nx*(y2-1))
                end

            end
        end

        table[c] = list

    end

    return table

end

############################################################
# Build spatial cell list
############################################################

function build_cell_grid!(sim::Simulation)

    agents = sim.agents
    grid = sim.grid

    fill!(grid.cell_count,0)
    fill!(grid.cell_offset,0)

    N = length(agents.x)

    ########################################
    # Pass 1: count agents per cell
    ########################################

    @threads for i in 1:N

        c = cell_index(
            agents.x[i],
            agents.y[i],
            sim.cell_size,
            sim.nx
        )

        @atomic grid.cell_count[c] += 1

    end

    ########################################
    # Prefix sum
    ########################################

    s = 1

    for c in 1:sim.ncells
        grid.cell_start[c] = s
        s += grid.cell_count[c]
    end

    ########################################
    # Pass 2: fill agent list
    ########################################

    @threads for i in 1:N

        c = cell_index(
            agents.x[i],
            agents.y[i],
            sim.cell_size,
            sim.nx
        )

        idx = atomic_add!(grid.cell_offset,c,1)

        grid.agent_index[
            grid.cell_start[c] + idx
        ] = i

    end

end

############################################################
# Ecological interaction kernel
############################################################

@inline function interaction_kernel!(i,j,sim)

    agents = sim.agents

    ti = agents.trait[i]
    tj = agents.trait[j]

    d = ti - tj
    ad = abs(d)

    ########################################
    # Competition
    ########################################

    comp = exp(-(ad^2)/(2f0*sim.competition_sigma^2))

    agents.energy[i] -= 0.01f0 * comp
    agents.energy[j] -= 0.01f0 * comp

    ########################################
    # Mating
    ########################################

    mate = exp(-(ad^2)/(2f0*sim.mating_sigma^2))

    if mate > 0.8f0
        agents.energy[i] += 0.05f0
        agents.energy[j] += 0.05f0
    end

    ########################################
    # Predation
    ########################################

    if d > sim.predation_threshold
        agents.energy[i] += 0.1f0
        agents.energy[j] -= 0.1f0

    elseif d < -sim.predation_threshold
        agents.energy[j] += 0.1f0
        agents.energy[i] -= 0.1f0
    end

end

############################################################
# Interaction step
############################################################

function interaction_step!(sim::Simulation)

    agents = sim.agents
    grid = sim.grid

    N = length(agents.x)

    @threads for i in 1:N

        xi = agents.x[i]
        yi = agents.y[i]

        c = cell_index(xi,yi,sim.cell_size,sim.nx)

        for nc in sim.neighbor_cells[c]

            start = grid.cell_start[nc]
            stop = start + grid.cell_count[nc] - 1

            for k in start:stop

                j = grid.agent_index[k]

                if j ≤ i
                    continue
                end

                dx = agents.x[j] - xi
                dy = agents.y[j] - yi

                r2 = dx*dx + dy*dy

                if r2 < sim.interaction_radius2
                    interaction_kernel!(i,j,sim)
                end

            end

        end

    end

end

############################################################
# Movement step
############################################################

function movement_step!(sim::Simulation)

    agents = sim.agents
    world = sim.world_size

    N = length(agents.x)

    @threads for i in 1:N

        agents.x[i] += 0.1f0*(rand(Float32)-0.5f0)
        agents.y[i] += 0.1f0*(rand(Float32)-0.5f0)

        agents.x[i] = mod(agents.x[i],world)
        agents.y[i] = mod(agents.y[i],world)

    end

end

############################################################
# One simulation step
############################################################

function step!(sim::Simulation)

    build_cell_grid!(sim)

    interaction_step!(sim)

    movement_step!(sim)

end

############################################################
# Initialization
############################################################

function init_simulation(N,world_size,cell_size)

    nx = Int(world_size/cell_size)
    ny = nx
    ncells = nx*ny

    agents = Agents(

        rand(Float32,N).*world_size,
        rand(Float32,N).*world_size,

        rand(Float32,N),
        ones(Float32,N)

    )

    neighbor_table = build_neighbor_table(nx,ny)

    grid = CellGrid(
        zeros(Int,ncells),
        zeros(Int,ncells),
        zeros(Int,ncells),
        Vector{Int}(undef,N)
    )

    return Simulation(

        agents,

        world_size,
        cell_size,

        nx,
        ny,
        ncells,

        neighbor_table,

        grid,

        1f0,    # interaction radius^2
        0.2f0,  # competition width
        0.05f0, # mating width
        0.3f0   # predation threshold

    )

end

