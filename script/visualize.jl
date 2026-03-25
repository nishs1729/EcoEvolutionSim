using Revise
using EcoEvolutionSim
using Random
using GLMakie

# 1. Initialize Simulation
sim = init_simulation("config.toml")

# 2. Setup Observables
points = Observable(Point2f.(sim.agents.x, sim.agents.y))
energies = Observable(sim.agents.energy .* 10.0f0) # Scale energy for markersize

color_trait_name = "fecundity"
color_trait = Observable(sim.agents.traits[color_trait_name])

# 3. Setup Figure
fig = Figure(size = (800, 800))
ax = Axis(fig[1, 1], title = "Eco-Evolutionary Simulation (Reflective Boundaries)",
          aspect = DataAspect(), 
          limits = (0, sim.config.world_size, 0, sim.config.world_size))

# 4. Plotting
# Agents
scatter!(ax, points,
         color = color_trait, 
         colormap = :viridis, 
         colorrange = (TRAIT_SPECS[color_trait_name].mean - 3.0f0 * TRAIT_SPECS[color_trait_name].sigma, 
                       TRAIT_SPECS[color_trait_name].mean + 3.0f0 * TRAIT_SPECS[color_trait_name].sigma),
         markersize = energies)

# Colorbar for traits
Colorbar(fig[1, 2], label = "Trait Value", colormap = :viridis, limits = (0, 1))

display(fig)

# 6. Update function
function update_visuals!(sim, points, color_trait, energies)
    points[] = Point2f.(sim.agents.x, sim.agents.y)
    color_trait[] = sim.agents.traits[color_trait_name]
    energies[] = sim.agents.energy .* 10.0f0
end

# 7. Animation Loop
try
    while isopen(fig.scene)
        for _ in 1:sim.config.steps_per_frame
            step!(sim)
        end
        update_visuals!(sim, points, color_trait, energies)
        sleep(0.01)
    end
catch e
    if !(e isa InterruptException)
        rethrow(e)
    end
end

