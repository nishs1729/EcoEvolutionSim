using Revise
using EcoEvolutionSim
using Random
using GLMakie

# 1. Initialize Simulation
traits_path = "script/traits.toml"
sim = init_simulation("config.toml", traits_path)
trait_specs = load_traits(traits_path)

# 2. Setup Observables
points = Observable(Point2f.(sim.agents.x, sim.agents.y))
energies = Observable(sim.agents.energy .* 10.0f0) # Scale energy for markersize
age = Observable(sim.agents.age) # Scale energy for markersize

color_trait_name = "fecundity"
color_trait = Observable(sim.agents.traits[Symbol(color_trait_name)])

# 3. Setup Figure
fig = Figure(size = (800, 800))
ax = Axis(fig[1, 1], title = "Eco-Evolutionary Simulation (Reflective Boundaries)",
          aspect = DataAspect(), 
          limits = (0, sim.config.world_size, 0, sim.config.world_size))

# 4. Plotting
# Agents
spec = trait_specs[color_trait_name]
scatter!(ax, points,
         color = color_trait, 
         colormap = :viridis, 
         colorrange = (spec.mean - 3.0f0 * spec.sigma, 
                       spec.mean + 3.0f0 * spec.sigma),
         markersize = age)
        #  markersize = energies)

# Colorbar for traits
Colorbar(fig[1, 2], label = "Trait Value", colormap = :viridis, limits = (0, 1))

display(fig)

# 6. Update function
function update_visuals!(sim, points, color_trait, energies)
    alive = sim.agents.alive

    points[] = Point2f.(sim.agents.x[alive], sim.agents.y[alive])
    color_trait[] = sim.agents.traits[Symbol(color_trait_name)][alive]
    # energies[] = sim.agents.energy .* 10.0f0
    age[] = sim.agents.age[alive]
end

# 7. Animation Loop
try
    while isopen(fig.scene)
        for _ in 1:sim.config.steps_per_frame
            step!(sim)
            # print(sim.agents.age[1], " ", sim.agents.energy[1], "\n")
        end
        update_visuals!(sim, points, color_trait, energies)
        sleep(0.001)
    end
catch e
    if !(e isa InterruptException)
        rethrow(e)
    end
end
