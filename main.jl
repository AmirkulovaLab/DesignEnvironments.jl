using DesignEnvironments
using ReinforcementLearning
ENV["GKSwstype"] = "nul"

env = CylinderEnv(M=10, continuous=true, grid_size=10.0)
policy = RandomPolicy(action_space(env))

render(env, policy, path="thing.mp4")
