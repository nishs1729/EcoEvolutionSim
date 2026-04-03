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

# Spatial Iterator
export nearby_agents

struct NearbyAgentIterator{S}
    sim::S
    agent_id::Int
    r_sq::Float32
end

function nearby_agents(sim::Simulation, agent_id::Int, r::Float32)
    NearbyAgentIterator{typeof(sim)}(sim, agent_id, r * r)
end

Base.IteratorSize(::Type{<:NearbyAgentIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:NearbyAgentIterator}) = Base.HasEltype()
Base.eltype(::Type{<:NearbyAgentIterator}) = Int

function Base.iterate(iter::NearbyAgentIterator)
    sim = iter.sim
    i = iter.agent_id
    env = sim.env
    x, y = sim.agents.x[i], sim.agents.y[i]
    c = cell_index(x, y, env.inv_cell_size, env.nx, env.ny)
    return iterate_next(iter, 1, 1, c)
end

function Base.iterate(iter::NearbyAgentIterator, state::Tuple{Int, Int, Int})
    return iterate_next(iter, state[1], state[2], state[3])
end

@inline function iterate_next(iter::NearbyAgentIterator, n_idx::Int, j_idx::Int, origin_c::Int)
    sim = iter.sim
    i = iter.agent_id
    env = sim.env
    r_sq = iter.r_sq
    x0, y0 = sim.agents.x[i], sim.agents.y[i]
    
    @inbounds while n_idx <= env.neighbors.count[origin_c]
        c = env.neighbors.neighbors[n_idx, origin_c]
        count = env.grid.cell_count[c]
        start = env.grid.cell_start[c]
        
        while j_idx <= count
            j = env.grid.agent_index[start + j_idx - 1]
            j_idx += 1
            if j != i && sim.agents.alive[j]
                dx = sim.agents.x[j] - x0
                dy = sim.agents.y[j] - y0
                if dx*dx + dy*dy <= r_sq
                    return (j, (n_idx, j_idx, origin_c))
                end
            end
        end
        n_idx += 1
        j_idx = 1
    end
    return nothing
end