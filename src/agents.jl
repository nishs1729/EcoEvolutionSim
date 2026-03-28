struct Agents
    x::Vector{Float32}
    y::Vector{Float32}
    vx::Vector{Float32}
    vy::Vector{Float32}
    gender::Vector{UInt8} # 0: female, 1: male
    energy::Vector{Float32}
    age::Vector{Float32}
    alive::Vector{Bool}
    traits::Dict{String, Vector{Float32}}
end

function Agents(N, world_size, traits_dict::Dict{String, Vector{Float32}})
    Agents(
        rand(Float32, N) .* world_size,
        rand(Float32, N) .* world_size,
        zeros(Float32, N),
        zeros(Float32, N),
        rand(UInt8(0):UInt8(1), N),
        rand(Float32, N),
        rand(Float32, N) .* 10f0,   # initial age
        trues(N),                   # alive
        traits_dict
    )
end

function age_step!(sim::Simulation)
    agents = sim.agents
    age = agents.age
    alive = agents.alive

    Threads.@threads for i in eachindex(age)
        if alive[i]
            age[i] += 0.01f0
        end
    end
end

function death_step!(sim::Simulation)
    agents = sim.agents
    age = agents.age
    alive = agents.alive

    μ0 = 0.0001f0
    k = 0.05f0

    Threads.@threads for i in eachindex(age)
        if !alive[i]
            continue
        end

        a = age[i]
        mortality = μ0 * exp(k*a)
        if rand(Float32) < mortality
            alive[i] = false
        end
    end
end