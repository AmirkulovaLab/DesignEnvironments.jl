export DesignEnvironment

abstract type AbstractDesign end
abstract type AbstractObjective end

const IS_CONTINUOUS = true
const EPISODE_LENGTH = 100
const PENALTY_WEIGHT = 0.1

"""
Base environment for episodic design optimization
"""
mutable struct DesignEnvironment <: AbstractEnv
    design::AbstractDesign
    objective::AbstractObjective
    is_continuous::Bool

    episode_length::Int
    penalty_weight::Float64

    timestep::Int
    penalty::Float64
end

function DesignEnvironment(;
    design::AbstractDesign,
    objective::AbstractObjective,
    is_continuous::Bool = IS_CONTINUOUS,
    episode_length::Int = EPISODE_LENGTH,
    penalty_weight::Float64 = PENALTY_WEIGHT)

    ## DesignEnvironment starting at timestep 0
    env = DesignEnvironment(
        design, objective, is_continuous, episode_length, penalty_weight, 0, 0.0
        )

    reset!(env)
    return env
end

function RLBase.reset!(env::DesignEnvironment)
    env.timestep = 0
    reset_design!(env.design)
    env.objective(env.design)
end

function RLBase.is_terminated(env::DesignEnvironment)
    return env.timestep >= env.episode_length
end

function RLBase.state_space(env::DesignEnvironment)
    return state_space(env.design, env.objective)
end

function RLBase.action_space(env::DesignEnvironment)
    return action_space(env.design, env.is_continuous)
end

function RLBase.state(env::DesignEnvironment)
    return vcat(now(env.design), now(env.objective))
end

function RLBase.reward(env::DesignEnvironment)
    return - (metric(env.objective) + env.penalty * env.penalty_weight)
end

function (env::DesignEnvironment)(action)
    env.timestep += 1
    env.penalty = 0.0
    ## make adjustments to design and record penalty
    env.penalty = env.design(action)
    env.objective(env.design)
end

function img(env::DesignEnvironment)
    return plot(img(env.design), img(env.objective), layout=@layout([a{0.6w} b]))
end