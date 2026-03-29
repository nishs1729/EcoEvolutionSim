include("dstructs.jl")

const MIN_PARALLEL_N = 1024

function Configurations.from_dict(::Type{MovementStrategy}, x::String)
    s = Symbol(uppercase(x))
    for i in instances(MovementStrategy)
        if Symbol(i) == s
            return i
        end
    end
    error("'$x' is not a valid MovementStrategy. Options: $(instances(MovementStrategy))")
end

Base.convert(::Type{MovementStrategy}, x::String) = Configurations.from_dict(MovementStrategy, x)

function load_config(path::String)
    return from_toml(Config, path)
end