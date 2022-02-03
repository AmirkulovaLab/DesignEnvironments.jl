function RLBase.action_space(design::Configuration, objective::TSCS, is_continuous::Bool)
    if is_continuous
        ## in the case of continuous actions the actino space will be a vector
        step_size = design.max_vel
        action_space = Space([-step_size..step_size for _ in 1:(2 * design.M)])
    else
        ## in the case of discrete actions it will be an integer within a range.
        action_space = Base.OneTo(4 * design.M)
    end

    return action_space
end

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