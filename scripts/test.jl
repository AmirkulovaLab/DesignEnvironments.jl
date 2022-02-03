using DesignEnvironments
using ReinforcementLearning
using Plots

env = DesignEnvironment(
    state_type = SequenceVectorState,
    is_continuous = false,
    episode_length = 100,
    penalty_weight = 0.1,
    design = Configuration(
        M = 10, 
        plane_size = 20.0,
        max_vel = 0.1,
        vel_decay = 0.7,
        min_distance = 0.1),
    objective = TSCS(
        k0amax = 1.5,
        k0amin = 0.5,
        nfreq = 15,
        rho = 1000.0,
        c0 = 1480.0,
        aa = 1.0))

# traj = CircularArraySARTTrajectory(
#     capacity = 1000, 
#     state = VectorState => ())

# policy = RandomPolicy(action_space(env))
# reset!(env)
# while !is_terminated(env)
#     a = policy(env)
#     push!(traj, state = state(env), action = a)
#     env(a)
#     push!(traj, state = state(env), reward = reward(env), terminal = is_terminated(env))
#     display(env.objective.Q_RMS)
# end

# savefig(plot(env), "plot")



