ENV["GKSwstype"] = "nul"

using DesignEnvironments
using ReinforcementLearning

env = CylinderEnv(M=10, grid_size=10.0, continuous=false, max_velocity=0.1)
policy = RandomPolicy(action_space(env))

DE.render(env, policy, path="new.mp4", max_tscs=10.0)
