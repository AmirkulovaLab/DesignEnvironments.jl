function RLBase.state_space(design::Configuration, objective::TSCS, ::Type{VectorState})
    M = design.M
    nfreq = objective.nfreq
    n_coords = M * 2
    
    return Space([-Inf..Inf for _ in 1:(3 * n_coords + nfreq + 1)])
end

"""
Obtains the current state of the configuration in the form of a `Vector`. The state
is comprised of the position and velocity of the cylinders.
"""
function RLBase.state(design::Configuration, objective::TSCS, ::Type{VectorState})
    pos = coords_to_x(design.pos)
    vel = coords_to_x(design.vel)
    Q = objective.Q
    qV = objective.qV
    Q_RMS = objective.Q_RMS

    features = vcat(pos, vel, Q, qV, Q_RMS)
    return VectorState(features)
end