function CellGrid(ncells, N)
    return CellGrid(
        zeros(Int, ncells),
        zeros(Int, ncells),
        zeros(Int, ncells),
        Vector{Int32}(undef, N)
    )
end

@inline function cell_index(x, y, cell_size, nx, ny)
    ix = clamp(Int(floor(x / cell_size)) + 1, 1, nx)
    iy = clamp(Int(floor(y / cell_size)) + 1, 1, ny)
    return ix + nx * (iy - 1)
end

############################################################
# Build spatial cell list
############################################################
function build_cell_grid!(sim::Simulation)
    agents = sim.agents
    env = sim.env
    ncells = env.ncells
    cell_size = env.cell_size
    nx = env.nx
    alive = agents.alive

    fill!(env.grid.cell_count, 0)
    fill!(env.grid.cell_offset, 0)

    N = length(agents.x)

    # For small N, serial is much faster than parallel with atomic contention
    @inbounds for i in 1:N
        if !alive[i]
            continue
        end
        c = cell_index(agents.x[i], agents.y[i], cell_size, nx, env.ny)
        env.grid.cell_count[c] += 1
    end

    s = 1
    @inbounds for c in 1:ncells
        env.grid.cell_start[c] = s
        s += env.grid.cell_count[c]
    end

    @inbounds for i in 1:N
        if !alive[i]
            continue
        end
        c = cell_index(agents.x[i], agents.y[i], cell_size, nx, env.ny)
        idx = env.grid.cell_offset[c]
        env.grid.cell_offset[c] += 1
        env.grid.agent_index[env.grid.cell_start[c] + idx] = i
    end
end

function build_neighbor_table(nx, ny)
    ncells = nx * ny
    table = Vector{Vector{Int}}(undef, ncells)
    for c in 1:ncells
        iy = (c - 1) ÷ nx + 1
        ix = c - nx * (iy - 1)
        list = Int[]
        for dy in -1:1
            for dx in -1:1
                x2 = ix + dx
                y2 = iy + dy
                if 1 ≤ x2 ≤ nx && 1 ≤ y2 ≤ ny
                    push!(list, x2 + nx * (y2 - 1))
                end
            end
        end
        table[c] = list
    end
    return table
end