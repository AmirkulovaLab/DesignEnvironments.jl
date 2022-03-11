using DesignEnvironments
using ReinforcementLearning
include("../src/core/multi_agent_design_environment.jl")

## define design
cylinders = Configuration(
    M = 4,
    plane_size = 10.0,
    max_vel = 0.05,
    vel_decay = 1.0,
    min_distance = 0.1)

tscs = TSCS(
    k0amax = 0.45,
    k0amin = 0.35,
    nfreq = 10,
    aa = maximum(cylinders.radii),
    rho = DE.RHO,
    c0 = DE.C0)

state_type = VectorState
## define environment
core = DesignEnvironment(
    design = cylinders,
    # objective = pressure_amplitude,
    objective = tscs,
    state_type = state_type,
    is_continuous = true,
    episode_length = 100,
    penalty_weight = 0.1)

env = MultiAgentDesignEnvironment(5, core)