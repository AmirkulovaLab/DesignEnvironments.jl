using DesignEnvironments
using ReinforcementLearning
using Plots
using ProgressMeter

M = 7

env = DesignEnvironment(
    design = CoreConfiguration(M = M, plane_size = 15.0, vel_decay=1.0),
    # design = Configuration(M = M, plane = Square(15.0)),
    objective = TSCS(k0amax = 1.0, k0amin = 0.3, nfreq = 15))

policy = RandomPolicy(action_space(env))

DE.render(env, policy, "anim.mp4")