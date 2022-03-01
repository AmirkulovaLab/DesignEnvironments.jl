function RLBase.action_space(design::Configuration, ::AbstractObjective, is_continuous::Bool)
    if is_continuous
        ## in the case of continuous actions the actino space will be a vector
        step_size = design.max_vel
        action_space = Space([-step_size..step_size for _ in 1:(2 * design.M)])
    else
        ## in the case of discrete actions it will be an integer within a range.
        action_space = 1:(4 * design.M)
    end

    return action_space
end