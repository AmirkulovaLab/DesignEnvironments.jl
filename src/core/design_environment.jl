export DesignEnvironment, AbstractDesign, AbstractObjective, AbstractState, stack

abstract type AbstractDesign end
abstract type AbstractObjective end
abstract type AbstractState end

"""
Stacks one or many states into an array
"""
function stack(s::AbstractState...)
    return hcat(extract.(s)...)
end

function RLCore.send_to_device(
        D::Union{CuDevice, Val{:cpu}}, 
        s::Union{AbstractState, Vector{AbstractState}})

    send_to_device(D, unpack(s))
end

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
mutable struct DesignEnvironment{D <: AbstractDesign, O <: AbstractObjective, S <: AbstractState} <: AbstractEnv
    design::D
    objective::O
    state_type::Type{S}
    is_continuous::Bool

    episode_length::Int
    penalty_weight::Float64

    timestep::Int
    penalty::Float64
end

RLBase.InformationStyle(::DesignEnvironment) = PERFECT_INFORMATION
RLBase.ChanceStyle(::DesignEnvironment) = DETERMINISTIC
RLBase.RewardStyle(::DesignEnvironment) = STEP_REWARD

function DesignEnvironment(;
    design::AbstractDesign,
    objective::AbstractObjective,
    state_type::Type,
    is_continuous::Bool,
    episode_length::Int,
    penalty_weight::Float64)

    ## DesignEnvironment starting at timestep 0
    env = DesignEnvironment(
        design, objective, state_type,
        is_continuous, episode_length, penalty_weight, 0, 0.0
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

function RLBase.reward(env::DesignEnvironment)
    return - (metric(env.objective) + env.penalty * env.penalty_weight)
end

function RLBase.action_space(env::DesignEnvironment)
    return action_space(env.design, env.objective, env.is_continuous)
end

function RLBase.state_space(env::DesignEnvironment)
    return state_space(env.design, env.objective, env.state_type)
end

function RLBase.state(env::DesignEnvironment)
    return state(env.design, env.objective, env.state_type)
end

function (env::DesignEnvironment)(action)
    ## increment timestep
    env.timestep += 1
    ## reset penalty
    env.penalty = 0.0
    ## make adjustments to design and record penalty
    env.penalty = env.design(action)
    ## calculate objective on new design
    env.objective(env.design)
end

function Plots.plot(env::DesignEnvironment, objective_scale::Tuple)
    return plot(
        plot(env.design), 
        plot(env.objective, objective_scale), 
        layout=@layout([a{0.6w} b]))
end

function Plots.plot(env::DesignEnvironment)
    return plot(env, scale(env.objective))
end