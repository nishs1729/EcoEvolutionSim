function load_traits(path::String)
    raw = TOML.parsefile(path)
    traits = Dict{String,TraitSpec}()
    for (name, spec) in raw
        traits[name] = TraitSpec(Float32(spec["mean"]), Float32(spec["sigma"]))
    end
    return traits
end

function initialize_traits(specs::Dict{String,TraitSpec}, max_agents, initial_agents)
    data = Dict{Symbol,Vector{Float32}}()

    for (name, spec) in specs
        v = Vector{Float32}(undef, max_agents)
        randn!(v)                                   # fill with N(0,1) in-place
        @inbounds @simd for i in eachindex(v)
            v[i] = spec.mean + spec.sigma * v[i]   # transform in-place, no extra alloc
        end
        data[Symbol(name)] = v
    end

    return NamedTuple(data)
end
