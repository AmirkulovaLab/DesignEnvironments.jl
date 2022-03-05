function RLBase.state_space(design::Configuration, objective::PressureAmplitude, ::Type{SequenceState})
    M = design.M
    nfreq = objective.nfreq

    seq_space = Space(Matrix(hcat([[-Inf..Inf for _ in 1:(4 + nfreq + 1)] for _ in 1:M]...)))
    return seq_space
end

function RLBase.state(design::Configuration, objective::PressureAmplitude, ::Type{SequenceState})
    pos = design.pos
    vel = design.vel
    Q = objective.Q
    Q_RMS = metric(objective)
    seq = hcat(pos, vel) |> transpose |> Matrix

    seq_len = size(seq, 2)
    features = hcat(fill(Q, seq_len)...)

    return SequenceState(vcat(seq, features))
end