const _STRATEGY_MAP = Dict(string(s) => s for s in instances(MovementStrategy))

function Configurations.from_dict(::Type{MovementStrategy}, x::String)
    s = uppercase(x)
    haskey(_STRATEGY_MAP, s) || error("'$x' is not a valid MovementStrategy. Options: $(instances(MovementStrategy))")
    return _STRATEGY_MAP[s]
end

Base.convert(::Type{MovementStrategy}, x::String) = Configurations.from_dict(MovementStrategy, x)

function load_config(path::String)
    return from_toml(Config, path)
end