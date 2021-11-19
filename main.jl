ENV["GKSwstype"] = "nul"
using DesignEnvironments
using ReinforcementLearning
using Plots

env = CylinderEnv(
    M=4,
    grid_size=5.0,
    k0amax=1.0,
    continuous=false,
    step_size=0.1,
    nfreq=10,
    velocity_decay=0.9,
    physics=true)

# policy = RandomPolicy(action_space(env))

mutable struct DummyPolicy <: AbstractPolicy
    action
    count::Int

    function DummyPolicy(env::CylinderEnv)
        action = env |> action_space |> rand
        return new(action, 0)
    end
end

function (policy::DummyPolicy)(env::CylinderEnv)
    if policy.count >= 1
        try
            policy.action = zeros(size(action_space(env)))
        finally
            policy.action = env |> action_space |> rand
        end

    end

    policy.count += 1
    return policy.action
end

policy = DummyPolicy(env)

DE.render(env, policy, path="new.mp4", max_tscs=30.0, fps=15)
