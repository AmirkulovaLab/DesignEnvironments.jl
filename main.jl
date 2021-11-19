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
    nfreq=20,
    velocity_decay=1.0,
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
        policy.action = zeros(size(action_space(env)))
    end

    policy.count += 1
    return policy.action
end

policy = DummyPolicy(env)

DE.render(env, policy, path="new.mp4", max_tscs=30.0, fps=15)

# a = Animation()
#
# while !is_terminated(env)
#     env(policy(env))
#
#     frame(a, img(env))
#
#     display(env.penalty)
# end

## convert collections of images into gif
# gif(a, "new.mp4", fps=15)
# closeall()

# run(policy, env)
