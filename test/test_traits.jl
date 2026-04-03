@testset "Traits Generation" begin
    @testset "Trait Loading" begin
        temp_traits = tempname() * ".toml"
        open(temp_traits, "w") do io
            write(io, """
            [size]
            mean = 10.0
            sigma = 2.0
            
            [speed]
            mean = 1.5
            sigma = 0.5
            """)
        end
        
        specs = load_traits(temp_traits)
        @test haskey(specs, "size")
        @test haskey(specs, "speed")
        @test specs["size"].mean == 10.0f0
        @test specs["size"].sigma == 2.0f0
        @test specs["speed"].mean == 1.5f0
        @test specs["speed"].sigma == 0.5f0
        
        rm(temp_traits)
    end
    
    @testset "Trait Initialization" begin
        specs = Dict{String, TraitSpec}(
            "size" => TraitSpec(10.0f0, 2.0f0),
            "speed" => TraitSpec(1.5f0, 0.5f0)
        )
        
        N = 10000
        traits = EcoEvolutionSim.initialize_traits(specs, N, N)
        
        @test typeof(traits) <: NamedTuple
        @test haskey(traits, :size)
        @test haskey(traits, :speed)
        
        @test length(traits.size) == N
        @test length(traits.speed) == N
        
        # Test distribution properties
        @test isapprox(sum(traits.size) / N, 10.0f0, atol=0.1)
        @test isapprox(sum(traits.speed) / N, 1.5f0, atol=0.05)
    end
end