using Revise
using EcoEvolutionSim

greet()

############################################################
# Run simulation
############################################################

N = 1_000

world_size = 100f0
cell_size = 1f0

println("Initializing simulation...")

sim = init_simulation(
    N,
    world_size,
    cell_size
)

steps = 100

println("Running simulation...")

for t in 1:steps

    step!(sim)

    if t % 10 == 0
        println("step ",t)
    end

end
