export SequenceState

export SequenceState

struct SequenceState <: AbstractState
    sequence::Matrix{Float32}
    n_features::Int
    seq_len::Int
end

function SequenceState(sequence::Matrix)
    n_features = size(sequence, 1)
    seq_len = size(sequence, 2)

    return SequenceState(sequence, n_features, seq_len)
end

function Base.ndims(state::SequenceState)
    return 0
end

function extract(s::SequenceState)
    return unsqueeze(s.sequence, 2)
end