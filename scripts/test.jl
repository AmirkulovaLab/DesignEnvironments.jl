using DesignEnvironments
using ReinforcementLearning
using Plots
using Random
using Flux

env = DesignEnvironment(
    state_type = SequenceVectorState,
    is_continuous = false,
    episode_length = 100,
    penalty_weight = 0.1,
    design = Configuration(M = 5, plane_size = 20.0, max_vel = 0.1, vel_decay = 0.7, min_distance = 0.1),
    objective = TSCS(k0amax = 1.5, k0amin = 0.5, nfreq = 15, rho = 1000.0, c0 = 1480.0, aa = 1.0))

traj = CircularArraySARTTrajectory(
    capacity = 1000,
    state = SequenceVectorState => (),
    # action = Vector => env |> action_space |> size
    )

sampler = NStepBatchSampler{SARTS}(
    γ=Float32(0.9),
    n=1,
    stack_size=1,
    batch_size=16)

policy = RandomPolicy(action_space(env))
reset!(env)
while !is_terminated(env)

    a = policy(env)

    push!(
        traj,
        state = state(env), 
        action = a)

    env(a)

    push!(
        traj,
        state = state(env), 
        reward = reward(env), 
        terminal = is_terminated(env))

    display(env.objective.Q_RMS)
end

D = Val{:cpu}()
idx, batch = sample(traj, sampler)
s, a, r, t, s_ = batch

s = send_to_device.(D, stack(s...))
a = send_to_device(D, a)
r = send_to_device(D, r)
t = send_to_device(D, t)
s_ = send_to_device.(D, stack(s_...))

include("model.jl")
h_size = 64
n_actions = size(action_space(env))[1]
display(n_actions)

model = DuelingNetwork(
    base = lstm_block(env, h_size),
    val = Dense(2 * h_size, 1),
    adv = Dense(2 * h_size, n_actions))