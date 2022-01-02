function RLBase.state_space(config::Configuration, tscs::TSCS)
    return Space([-Inf..Inf for _ in Base.OneTo(6 * config.M + tscs.nfreq)])
end

function RLBase.action_space(config::Configuration, is_continuous::Bool)
    if is_continuous
        ## in the case of continuous actions the actino space will be a vector
        step_size = config.max_vel
        action_space = Space([-step_size..step_size for _ in 1:(2 * config.M)])
    else
        ## in the case of discrete actions it will be an integer within a range.
        action_space = Base.OneTo(4 * config.M)
    end

    return action_space
end

function RLBase.state_space(design::CoreConfiguration, tscs::DE.TSCS)
    return state_space(design.config, tscs)
end

function RLBase.action_space(design::CoreConfiguration, is_continuous::Bool)
    return action_space(design.config, is_continuous)
end