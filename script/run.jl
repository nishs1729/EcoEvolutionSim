using EcoEvolutionSim
using HDF5
using Random

config_file = "config.toml"
nsteps = 2000
save_every = 100
outfile = "data/simulation_output.h5"

# Initialize simulation
sim = init_simulation(config_file)
max_agents = length(sim.agents.x)
nsaves = nsteps ÷ save_every

println("Max agents: ", max_agents)
println("Snapshots to store: ", nsaves)

# Create HDF5 file
h5open(outfile, "w") do file
    file["world_size"] = sim.config.world_size
    file["save_every"] = save_every

    # Preallocated datasets
    x_dset = create_dataset(file, "x", Float32,
        (max_agents, nsaves), chunk=(max_agents, 1), compress=3)

    y_dset = create_dataset(file, "y", Float32,
        (max_agents, nsaves), chunk=(max_agents, 1), compress=3)

    energy_dset = create_dataset(file, "energy", Float32,
        (max_agents, nsaves), chunk=(max_agents, 1), compress=3)

    age_dset = create_dataset(file, "age", Float32,
        (max_agents, nsaves), chunk=(max_agents, 1), compress=3)

    pop_dset = create_dataset(file, "population", Int32,
        (nsaves,), chunk=(min(nsaves, 1000),))

    time_dset = create_dataset(file, "time", Int32,
        (nsaves,), chunk=(min(nsaves, 1000),))

    # Trait datasets
    trait_dsets = Dict{String, HDF5.Dataset}()

    for (name, _) in pairs(sim.agents.traits)
        sname = string(name)
        trait_dsets[sname] = create_dataset(
            file,
            "trait_" * sname,
            Float32,
            (max_agents, nsaves),
            chunk=(max_agents, 1),
            compress=3
        )
    end

    # Simulation loop
    snapshot = 0
    for step in 1:nsteps

        step!(sim)
        pop = sum(sim.agents.alive)

        if pop == 0
            println("Extinction reached at step $step. Stopping.")
            break
        end

        if step % save_every == 0

            snapshot += 1
            alive = sim.agents.alive

            # store data
            x_dset[1:pop, snapshot] = sim.agents.x[alive]
            y_dset[1:pop, snapshot] = sim.agents.y[alive]
            energy_dset[1:pop, snapshot] = sim.agents.energy[alive]
            age_dset[1:pop, snapshot] = sim.agents.age[alive]

            for (name, values) in pairs(sim.agents.traits)
                trait_dsets[string(name)][1:pop, snapshot] = values[alive]
            end

            pop_dset[snapshot] = pop
            time_dset[snapshot] = step

            println("Saved snapshot $snapshot (step $step, pop $pop)")
        end
    end

    println("\n=== Performance Report ===")
    show(to)
    println("\n\nSimulation complete.")
end