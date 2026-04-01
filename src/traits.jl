function load_traits(path::String)
    raw = TOML.parsefile(path)
    traits = Dict{String, TraitSpec}()
    for (name, spec) in raw
        traits[name] = TraitSpec(Float32(spec["mean"]), Float32(spec["sigma"]))
    end
    return traits
end

function initialize_traits(specs::Dict{String, TraitSpec}, n_agents)
    data = Dict{Symbol, Vector{Float32}}()
    buffer = Vector{Float32}(undef, n_agents)

    for (name, spec) in specs
        randn!(buffer)
        v = similar(buffer)
        @inbounds @simd for i in eachindex(buffer)
            v[i] = spec.mean + spec.sigma * buffer[i]
        end
        data[Symbol(name)] = v
    end

    return NamedTuple(data)
end
