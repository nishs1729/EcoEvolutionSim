module Traits

using ..TraitRegistry: register_trait
using Random

export continuous_trait_constructor, bounded_trait_constructor

function continuous_trait_constructor(N::Int)
    return rand(Float32, N)
end

function bounded_trait_constructor(N::Int, min_val::Float32, max_val::Float32)
    return rand(Float32, N) .* (max_val - min_val) .+ min_val
end

function __init__()
    register_trait("continuous", continuous_trait_constructor)
    register_trait("bounded", bounded_trait_constructor)
end

end
