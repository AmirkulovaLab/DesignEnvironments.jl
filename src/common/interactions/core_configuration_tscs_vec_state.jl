function RLBase.action_space(design::CoreConfiguration, objective::TSCS, is_continuous::Bool)
    return nothing
end

function RLBase.state_space(design::CoreConfiguration, objective::TSCS, state_type::Type{VectorState})
    return nothing
end

function RLBase.state(design::CoreConfiguration, objective::TSCS, state_type::Type{VectorState})
    return nothing
end

# function RLBase.state_space(design::CoreConfiguration, tscs::TSCS)
    
#     return state_space(merge_configs(design.core, design.config), tscs)
# end

# function RLBase.action_space(design::CoreConfiguration, is_continuous::Bool)
#     return action_space(design.config, is_continuous)
# end