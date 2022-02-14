export SequenceVectorState, stack

struct SequenceVectorState <: AbstractState
    seq_state::SequenceState
    vec_state::VectorState
end

function SequenceVectorState(seq::Matrix, features::Vector)
    return SequenceVectorState(SequenceState(seq), VectorState(features))
end

function stack(s::SequenceVectorState...)
    seq_state = stack(getfield.(s, :seq_state)...)
    vec_state = stack(getfield.(s, :vec_state)...)
    return (seq_state, vec_state)
end

# function RLCore.send_to_device(D::Val{:cpu}, s::AbstractArray{SequenceVectorState})
#     display("DIFOPJSJSF")
#     seq, features = stack(s...)
#     return (send_to_device(D, seq), send_to_device(D, features))
# end