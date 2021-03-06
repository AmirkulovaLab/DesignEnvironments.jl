function RLBase.state_space(design::Configuration, objective::TSCS, ::Type{SequenceState})
    M = design.M
    nfreq = objective.nfreq

    feature_dim = 6 + nfreq + 1

    return Matrix(hcat([[-Inf..Inf for _ in 1:feature_dim] for _ in 1:M]...))
end

"""
"""
function RLBase.state(design::Configuration, objective::TSCS, ::Type{SequenceState})

    scattering = vcat(objective.Q, DE.metric(objective))
    scattering = repeat(scattering, inner = (1, design.M))

    coords = hcat(
        design.pos, 
        design.vel, 
        DE.x_to_coords(objective.qV)) |> transpose |> Matrix

    return SequenceState(vcat(coords, scattering))
end