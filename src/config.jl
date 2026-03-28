@enum MovementStrategy begin
    RANDOM_WALK = 1
    LANGEVIN = 2
    CORRELATED_RW = 3
    ACTIVE_BROWNIAN = 4
end

@option struct Config
    # Simulation
    n_agents::Int = 500
    world_size::Float32 = 50.0f0
    cell_size::Float32 = 5.0f0
    steps_per_frame::Int = 10

    # Movement
    strategy::MovementStrategy = RANDOM_WALK
    base_speed::Float32 = 0.1f0

    # Ecology / interactions
    r_interact::Float32 = 1.0f0
end

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