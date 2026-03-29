## User defined trait specifications
const TRAIT_SPECS = Dict(
    "dispersal" => TraitSpec(1.0f0, 0.2f0),
    "fecundity" => TraitSpec(5.0f0, 1.0f0)
)

function initialize_traits(specs::Dict{String, TraitSpec}, n_agents)
    traits = Dict{String, Vector{Float32}}()
    buffer = Vector{Float32}(undef, n_agents)

    for (name, spec) in specs
        randn!(buffer)
        v = similar(buffer)
        @inbounds @simd for i in eachindex(buffer)
            v[i] = spec.mean + spec.sigma*buffer[i]
        end
        traits[name] = v
    end

    traits
end