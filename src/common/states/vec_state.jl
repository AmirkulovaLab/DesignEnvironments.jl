export VectorState, stack

"""
A design state which is represented by a vector.
"""
struct VectorState <: AbstractState
    features::Vector{Float32}
end

"""
Extract the vector information held within the state and add a dimention for stacking.
"""
function extract(s::VectorState)
    return unsqueeze(s.features, 2)
end

"""
Stacks one or many states into an array
"""
function stack(s::VectorState...)
    hcat(extract.(s)...)
end