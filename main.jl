ENV["GKSwstype"] = "nul"
using DesignEnvironments
using ReinforcementLearning
using Plots

env = CylinderEnv(
    M=10,
    grid_size=8.0,
    k0amax=1.0,
    continuous=true,
    step_size=0.1,
    nfreq=2,
    velocity_decay=1.0,
    physics=false)

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
        if env.continuous
            policy.action = zeros(size(action_space(env)))
        else
            policy.action = 0
        end
    end

    policy.count += 1
    return policy.action
end

policy = DummyPolicy(env)

DE.render(env, policy, path="new.mp4", max_tscs=30.0, fps=15)
