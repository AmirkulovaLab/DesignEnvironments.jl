export SequenceVectorState

struct SequenceVectorState <: AbstractState
    seq_state::SequenceState
    vec_state::VectorState
end

function SequenceVectorState(seq::Matrix, features::Vector)
    return SequenceVectorState(SequenceState(seq), VectorState(features))
end

function DE.stack(s::SequenceVectorState...)
    seq_state = DE.stack(getfield.(s, :seq_state)...)
    vec_state = DE.stack(getfield.(s, :vec_state)...)
    return (seq_state, vec_state)
end