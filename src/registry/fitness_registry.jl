module FitnessRegistry

export FITNESS_REGISTRY, register_fitness

const FITNESS_REGISTRY = Dict{String, Function}()

function register_fitness(name::String, constructor::Function)
    FITNESS_REGISTRY[name] = constructor
end

end
