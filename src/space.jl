function CellGrid(ncells, N)
    return CellGrid(
        zeros(Int, ncells),
        zeros(Int, ncells),
        zeros(Int, ncells),
        Vector{Int32}(undef, N)
    )
end

@inline function cell_index(x, y, inv_cell_size, nx, ny)
    ix = clamp(Int(floor(x * inv_cell_size)) + 1, 1, nx)
    iy = clamp(Int(floor(y * inv_cell_size)) + 1, 1, ny)
    return ix + nx * (iy - 1)
end

############################################################
# Build spatial cell list
############################################################
function build_cell_grid!(sim::Simulation)
    agents = sim.agents
    env = sim.env
    ncells = env.ncells
    inv_cell_size = env.inv_cell_size
    nx = env.nx
    alive = agents.alive

    fill!(env.grid.cell_count, 0)
    fill!(env.grid.cell_offset, 0)

    N = agents.max_id
    max_tid = Threads.maxthreadid()
    # max_tid = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads() # version-safe check for Threads.maxthreadid, ensuring compatibility with older Julia versions while still correctly handling thread-local arrays in more recent versions.

    if Threads.nthreads() == 1
        @inbounds for i in 1:N
            if !alive[i]
                continue
            end
            c = cell_index(agents.x[i], agents.y[i], inv_cell_size, nx, env.ny)
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
            c = cell_index(agents.x[i], agents.y[i], inv_cell_size, nx, env.ny)
            idx = env.grid.cell_offset[c]
            env.grid.cell_offset[c] += 1
            env.grid.agent_index[env.grid.cell_start[c] + idx] = i
        end
        return
    end

    # Parallel version using thread-local accumulation to avoid atomic contention
    counts = zeros(Int32, ncells, max_tid)
    
    Threads.@threads :static for i in 1:N
        if !alive[i]
            continue
        end
        c = cell_index(agents.x[i], agents.y[i], inv_cell_size, nx, env.ny)
        @inbounds counts[c, Threads.threadid()] += 1
    end

    Threads.@threads for c in 1:ncells
        s = Int32(0)
        for t in 1:max_tid
            @inbounds s += counts[c, t]
        end
        @inbounds env.grid.cell_count[c] = s
    end

    s = Int32(1)
    @inbounds for c in 1:ncells
        env.grid.cell_start[c] = s
        s += env.grid.cell_count[c]
    end

    offsets = zeros(Int32, ncells, max_tid)
    Threads.@threads for c in 1:ncells
        @inbounds offsets[c, 1] = 0
        for t in 2:max_tid
            @inbounds offsets[c, t] = offsets[c, t-1] + counts[c, t-1]
        end
    end

    Threads.@threads :static for i in 1:N
        if !alive[i]
            continue
        end
        c = cell_index(agents.x[i], agents.y[i], inv_cell_size, nx, env.ny)
        tid = Threads.threadid()
        @inbounds begin
            idx = offsets[c, tid]
            offsets[c, tid] += 1
            env.grid.agent_index[env.grid.cell_start[c] + idx] = i
        end
    end
end

function build_neighbor_table(nx, ny)
    ncells = nx * ny
    neighbors = Matrix{Int32}(undef, 9, ncells)
    count = zeros(Int8, ncells)
    
    Threads.@threads for c in 1:ncells
        iy = (c - 1) ÷ nx + 1
        ix = c - nx * (iy - 1)
        n_count = 0
        for dy in -1:1
            for dx in -1:1
                x2 = ix + dx
                y2 = iy + dy
                if 1 ≤ x2 ≤ nx && 1 ≤ y2 ≤ ny
                    n_count += 1
                    @inbounds neighbors[n_count, c] = x2 + nx * (y2 - 1)
                end
            end
        end
        @inbounds count[c] = n_count
    end
    return NeighborTable(neighbors, count)
end