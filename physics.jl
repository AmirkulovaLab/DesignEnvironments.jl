using DesignEnvironments
using ReinforcementLearning
ENV["GKSwstype"] = "nul"


env = Physics(CylinderEnv(M=20, grid_size=10.0, continuous=false))
# policy = RandomPolicy(action_space(env))
#
# while !is_terminated(env)
#     display(env.env.Q_RMS)
#     action = policy(env)
#     env(action)
# end

#
# render(env, policy, path="M=10.mp4")
#
# env.env.x = (2 * env.env.params.grid_size) .* (rand(Float64, 2 * env.env.params.M) .- 0.5)
# resolve_cylinder_collisions(env)

display(img(env.env))
