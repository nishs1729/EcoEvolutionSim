using EcoEvolutionSim
using BenchmarkTools
using Profile
using PProf
using Random

# Load a standard configuration
config_path = "config.toml"
traits_path = "script/traits.toml"

println("=== Benchmarking initialization ===")
@btime init_simulation($config_path, $traits_path)

sim = init_simulation(config_path, traits_path)

println("\n=== Benchmarking full step! ===")
@btime step!($sim)

println("\n=== Benchmarking build_cell_grid! ===")
@btime build_cell_grid!($sim)

println("\n=== Benchmarking ecology_step! ===")
@btime ecology_step!($sim)

println("\n=== Benchmarking movement_step! ===")
@btime movement_step!($sim)

println("\n=== Benchmarking trait access (NamedTuple) ===")
# Accessing first element of a trait
@btime $sim.agents.traits.fecundity[1]

# Summing a trait (vector operation)
@btime sum($sim.agents.traits.fecundity)

println("\n=== Starting Profile (1000 steps) ===")
Profile.clear()
@profile for i in 1:1000
    step!(sim)
end

println("Profile complete. Run `pprof()` in a REPL if needed, or exporting to file...")
# We can't easily open a browser here, so we'll just save the profile.
pprof(out="data/profile.pb.gz")
println("Profile saved to data/profile.pb.gz")
