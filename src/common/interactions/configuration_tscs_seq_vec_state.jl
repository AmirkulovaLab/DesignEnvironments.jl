function RLBase.state_space(design::Configuration, objective::TSCS, ::Type{SequenceVectorState})
    M = design.M
    nfreq = objective.nfreq

    seq_space = Space(Matrix(hcat([[-Inf..Inf for _ in 1:6] for _ in 1:M]...)))
    feature_space = Space([-Inf..Inf for _ in 1:(nfreq + 1)])

    return (seq_space, feature_space)
end


function RLBase.state(design::Configuration, objective::TSCS, ::Type{SequenceVectorState})
    pos = design.pos
    vel = design.vel
    qV = x_to_coords(objective.qV)
    Q = objective.Q
    Q_RMS = objective.Q_RMS

    seq = hcat(pos, vel, qV) |> Matrix |> transpose |> Matrix
    features = vcat(Q, Q_RMS)

    return SequenceVectorState(SequenceState(seq), VectorState(features))
end