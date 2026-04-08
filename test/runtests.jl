using EcoEvolutionSim
using Test
using Random
using TOML
using Distributions

# Load the interaction functions that live in script/ so test_interactions can use them
include(joinpath(@__DIR__, "..", "script", "interactions.jl"))

@testset "EcoEvolutionSim.jl" begin
    include("test_config.jl")
    include("test_traits.jl")
    include("test_space.jl")
    include("test_movement.jl")
    include("test_ecology.jl")
    include("test_interactions.jl")
end