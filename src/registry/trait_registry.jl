module TraitRegistry

export TRAIT_REGISTRY, register_trait

const TRAIT_REGISTRY = Dict{String, Function}()

function register_trait(name::String, constructor::Function)
    TRAIT_REGISTRY[name] = constructor
end

end
