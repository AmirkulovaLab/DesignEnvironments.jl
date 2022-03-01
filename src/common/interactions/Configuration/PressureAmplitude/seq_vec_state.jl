function RLBase.state_space(design::Configuration, objective::PressureAmplitude, ::Type{SequenceVectorState})
    M = design.M
    nfreq = objective.nfreq

    seq_space = Space(Matrix(hcat([[-Inf..Inf for _ in 1:4] for _ in 1:M]...)))
    feature_space = Space([-Inf..Inf for _ in 1:(nfreq + 1)])

    return (seq_space, feature_space)
end

function RLBase.state(design::Configuration, objective::PressureAmplitude, ::Type{SequenceVectorState})
    pos = design.pos
    vel = design.vel
    Q = objective.Q
    Q_RMS = metric(objective)
    seq = hcat(pos, vel) |> transpose |> Matrix
    features = vcat(Q, Q_RMS)

    return SequenceVectorState(SequenceState(seq), VectorState(features))
end