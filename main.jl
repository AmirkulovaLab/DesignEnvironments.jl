ENV["GKSwstype"] = "nul"
using DesignEnvironments
using ReinforcementLearning
using Plots

env = CylinderEnv(
    M=10,
    grid_size=8.0,
    k0amax=1.0,
    continuous=false,
    step_size=0.1,
    nfreq=30,
    physics=true)

policy = RandomPolicy(action_space(env))

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
