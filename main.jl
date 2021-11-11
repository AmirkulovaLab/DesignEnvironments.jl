ENV["GKSwstype"] = "nul"

using DesignEnvironments
using ReinforcementLearning

env = CylinderEnv(
    M=5,
    grid_size=5.0,
    k0amax=2.0,
    continuous=false,
    max_velocity=0.1,
    step_size=0.1,
    physics=true)

policy = RandomPolicy(action_space(env))

# run(policy, env)
DE.render(env, policy, path="new.mp4", max_tscs=20.0)
