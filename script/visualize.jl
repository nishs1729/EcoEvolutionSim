using Revise
using EcoEvolutionSim
using GLMakie
using LinearAlgebra

# 1. Initialize Simulation
config = load_config("config.toml")
sim = init_simulation(config)
world_size = config.world_size

# 2. Setup Observables
points = Observable(Point2f.(sim.agents.x, sim.agents.y))
traits = Observable(sim.agents.trait)
energies = Observable(sim.agents.energy .* 10.0f0) # Scale energy for markersize
interaction_segments = Observable(Point2f[])

# 3. Setup Figure
fig = Figure(size = (800, 800))
ax = Axis(fig[1, 1], title = "Eco-Evolutionary Simulation (Reflective Boundaries)",
          aspect = DataAspect(), 
          limits = (0, world_size, 0, world_size))

# 4. Plotting
# Interaction lines (like bonds)
linesegments!(ax, interaction_segments, color = (:gray, 0.2), linewidth = 1)

# Agents
scatter!(ax, points, 
         color = traits, 
         colormap = :viridis, 
         colorrange = (0, 1),
         markersize = 10)

# Colorbar for traits
Colorbar(fig[1, 2], label = "Trait Value", colormap = :viridis, limits = (0, 1))

display(fig)

# 5. Interaction detection for visualization
function get_interaction_segments(sim)
    segments = Point2f[]
    agents = sim.agents
    grid = sim.grid
    N = length(agents.x)
    
    for i in 1:N
        xi, yi = agents.x[i], agents.y[i]
        c = cell_index(xi, yi, sim.cell_size, sim.nx)
        
        # Check current and neighbor cells
        for nc in sim.neighbor_cells[c]
            start = grid.cell_start[nc]
            stop = start + grid.cell_count[nc][] - 1
            
            for k in start:stop
                j = grid.agent_index[k]
                if j <= i continue end
                
                dx = agents.x[j] - xi
                dy = agents.y[j] - yi
                
                r2 = dx*dx + dy*dy
                if r2 < sim.interaction_radius2
                    push!(segments, Point2f(agents.x[i], agents.y[i]))
                    push!(segments, Point2f(agents.x[j], agents.y[j]))
                end
            end
        end
    end
    return segments
end

# 6. Update function
function update_visuals!(sim, points, traits, energies, interaction_segments)
    points[] = Point2f.(sim.agents.x, sim.agents.y)
    traits[] = sim.agents.trait
    energies[] = sim.agents.energy
    interaction_segments[] = get_interaction_segments(sim)
end

# 7. Animation Loop
try
    while isopen(fig.scene)
        # Multiple steps per frame for smoother/faster simulation
        for _ in 1:2
            step!(sim)
        end
        update_visuals!(sim, points, traits, energies, interaction_segments)
        sleep(0.01)
    end
catch e
    if !(e isa InterruptException)
        rethrow(e)
    end
end
