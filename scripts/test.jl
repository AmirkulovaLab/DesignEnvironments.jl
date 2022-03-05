using DesignEnvironments
using ReinforcementLearning
using IntervalSets
using DesignEnvironments: metric
include("../src/common/interactions/Configuration/PressureAmplitude/seq_state.jl")

## define design
cylinders = Configuration(
    M = 8,
    plane_size = 10.0,
    max_vel = 0.7,
    vel_decay = 0.9,
    min_distance = 0.1)

## define objectives
pressure_amplitude = PressureAmplitude(
    xf = [12.0, 0.0],
    k0amax = 0.45,
    k0amin = 0.35,
    nfreq = 11,
    a = maximum(cylinders.radii),
    R2 = cylinders.plane.size,
    rho = DE.RHO,
    c0 = DE.C0,
    use_cuda = false)

## define environment
env = DesignEnvironment(
    design = cylinders,
    objective = pressure_amplitude,
    state_type = SequenceState,
    is_continuous = false,
    episode_length = 100,
    penalty_weight = 0.1)

s = state_space(env)