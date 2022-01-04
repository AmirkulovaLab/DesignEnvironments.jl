export DesignEnvironment

abstract type AbstractDesign end
abstract type AbstractObjective end

const IS_CONTINUOUS = true
const EPISODE_LENGTH = 100
const PENALTY_WEIGHT = 0.1

"""
Base environment for episodic design optimization.

# Example
```
env = DesignEnvironment(
    design = Configuration(M = M, plane = Square(10.0), vel_decay=0.9),
    objective = TSCS(k0amax = 1.0, k0amin = 0.3, nfreq = 15)
    )
```

# Parameters
- `design::AbstractDesign`
- `objective::AbstractObjective`
- `is_continuous::Bool`: specifies if design actions are continuous or discrete
- `episode_length::Int`: number of steps (actions) in a design episode
- `penalty_weight::Float64`: scalar which determines penalty for actions which result in invalid states
- `timestep::Int`: current step number in the episode
- `penalty::Float64`: holds the current penalty
"""
mutable struct DesignEnvironment{D <: AbstractDesign, O <: AbstractObjective} <: AbstractEnv
    design::D
    objective::O
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

function Plots.plot(env::DesignEnvironment, objective_scale::Tuple)
    return plot(
        plot(env.design), 
        plot(env.objective, objective_scale), 
        layout=@layout([a{0.6w} b]))
end

function Plots.plot(env::DesignEnvironment)
    return plot(
        plot(env.design), 
        plot(env.objective, scale(env.objective)),
        layout=@layout([a{0.6w} b]))
end