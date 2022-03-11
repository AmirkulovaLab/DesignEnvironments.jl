using DesignEnvironments
using ReinforcementLearning
using Plots

cylinders = Configuration(
    M = 8,
    plane_size = 10.0,
    max_vel = 0.2,
    vel_decay = 0.7,
    min_distance = 0.1
)

## define objectives
pressure_amplitude = PressureAmplitude(
    xf = [12.0, 0.0],
    k0amax = 2.0,
    k0amin = 0.35,
    nfreq = 100,
    a = maximum(cylinders.radii),
    R2 = cylinders.plane.size,
    rho = DE.RHO,
    c0 = DE.C0,
    use_cuda = false)

env = DesignEnvironment(
    design = cylinders,
    objective = pressure_amplitude,
    episode_length = 100,
    penalty_weight = 0.1,
    is_continuous = true,
    state_type = VectorState,
)

pressure_amplitude(cylinders)

savefig(plot(env), "pa")

