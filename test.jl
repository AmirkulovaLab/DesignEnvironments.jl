using DesignEnvironments
using ReinforcementLearning

env = CylinderEnv(M=10, grid_size=10.0, continuous=false)
policy = RandomPolicy(action_space(env))

DE.render(env, policy, path="new.mp4")
