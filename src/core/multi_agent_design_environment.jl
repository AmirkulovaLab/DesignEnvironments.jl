export MultiAgentDesignEnvironment

mutable struct MultiAgentDesignEnvironment <: AbstractEnv
    core::DesignEnvironment
    player::Int
end

function MultiAgentDesignEnvironment(core::DesignEnvironment)

end

RLBase.NumAgentStyle(::MultiAgentDesignEnvironment) = MultiAgent(2)
