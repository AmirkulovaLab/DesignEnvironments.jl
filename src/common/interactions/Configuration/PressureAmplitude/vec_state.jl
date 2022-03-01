function RLBase.state_space(design::Configuration, objective::PressureAmplitude, ::Type{VectorState})
    M = design.M
    nfreq = objective.nfreq
    n_coords = M * 2
    
    return Space([-Inf..Inf for _ in 1:(2 * n_coords + nfreq + 1)])
end

"""
Obtains the current state of the configuration in the form of a `Vector`. The state
is comprised of the position and velocity of the cylinders.
"""
function RLBase.state(design::Configuration, objective::PressureAmplitude, ::Type{VectorState})
    pos = coords_to_x(design.pos)
    vel = coords_to_x(design.vel)
    Q = objective.Q
    Q_RMS = metric(objective)

    features = vcat(pos, vel, Q, Q_RMS)
    return VectorState(features)
end