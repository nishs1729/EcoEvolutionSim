function CellGrid(ncells, N)
    n_threads = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads()
    return CellGrid(
        zeros(Int32, ncells),
        zeros(Int32, ncells),
        zeros(Int32, ncells),
        Vector{Int32}(undef, N),
        zeros(Int32, ncells, n_threads),
        zeros(Int32, ncells, n_threads),
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

    N = agents.max_id
    max_tid = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads()

    if Threads.nthreads() == 1
        fill!(env.grid.cell_count, Int32(0))
        fill!(env.grid.cell_offset, Int32(0))
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

    # Parallel version — 2 :static barriers, zero Task allocations
    #
    # Old design (4 barriers):
    #   :static count | dynamic reduce | serial prefix-sum | dynamic offsets | :static place
    #   The two dynamic @threads calls allocated Julia Tasks each step → ~77 KiB/call at 12t.
    #
    # New design (2 barriers):
    #   :static count | serial combined pass | :static place
    #   The serial pass does cell_start + cell_count + per-thread offsets in one sweep.
    #   Cost: O(ncells × max_tid) scalar ops — ~60k for typical configs, < 0.05 ms.
    counts  = env.grid.thread_counts
    offsets = env.grid.thread_offsets
    fill!(counts,  Int32(0))

    # Barrier 1: each thread counts agents into its private column of `counts`
    Threads.@threads :static for i in 1:N
        if !alive[i]
            continue
        end
        c = cell_index(agents.x[i], agents.y[i], inv_cell_size, nx, env.ny)
        @inbounds counts[c, Threads.threadid()] += 1
    end

    # Serial combined pass — replaces three separate phases:
    #   (a) reduce counts[:, 1:T] → cell_count[c]   (was a parallel @threads loop)
    #   (b) exclusive prefix-sum → cell_start[c]     (was a serial loop)
    #   (c) per-thread offsets within each cell       (was a parallel @threads loop)
    # All three are now one O(ncells × max_tid) pass with sequential memory access.
    s = Int32(1)
    @inbounds for c in 1:ncells
        env.grid.cell_start[c] = s
        running = Int32(0)
        for t in 1:max_tid
            offsets[c, t] = running          # exclusive prefix within this cell
            running += counts[c, t]
        end
        env.grid.cell_count[c] = running    # total agents in cell c
        s += running
    end

    # Barrier 2: each thread places its agents using the pre-computed offsets
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