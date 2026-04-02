using EcoEvolutionSim
using Test
using Random
using TOML

@testset "EcoEvolutionSim.jl" begin
    include("test_config.jl")
    include("test_traits.jl")
    include("test_space.jl")
    include("test_movement.jl")
    include("test_ecology.jl")
end