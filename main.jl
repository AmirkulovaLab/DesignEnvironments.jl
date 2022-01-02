using DesignEnvironments
using ReinforcementLearning
using Plots
using ProgressMeter

M = 7

env = DesignEnvironment(
    design = CoreConfiguration(M = M, plane_size = 15.0, vel_decay=0.9),
    # design = Configuration(M = M, plane = Square(15.0), vel_decay=0.8),
    objective = TSCS(k0amax = 1.0, k0amin = 0.3, nfreq = 15),
    is_continuous=false)

policy = RandomPolicy(action_space(env))

run(policy, env, StopWhenDone(), Render("anim.mp4"))