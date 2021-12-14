ENV["GKSwstype"] = "nul"
using DesignEnvironments
using Plots
using ReinforcementLearning

make_env() = CylinderEnv(
    M = 30,
    plane = Square(20.0),
    max_vel = 0.5,
    vel_decay = 1.0)

# envs = [make_env() for i in 1:10]

# policy = RandomPolicy(action_space(envs[1]))

MultiThreadEnv(make_env, 10)

# action = [policy(env) for env in envs]
# action[1]

# @time Threads.@threads for i in 1:10

#     reset!(env)
#     while !is_terminated(env)
#         env(policy(env))
#     end 
# end
