using DesignEnvironments
using ReinforcementLearning
using Plots

env = DesignEnvironment(

    ## environment arguments
    state_type = SequenceVectorState,
    is_continuous = false,
    episode_length = 100,
    penalty_weight = 0.1,

    ## constructing the design
    design = Configuration(
        M = 20,
        plane_size = 15.0,
        max_vel = 0.2,
        vel_decay = 0.8,
        min_distance = 0.1),

    ## constructing the objective
    objective = TSCS(
        k0amax = 1.0,
        k0amin = 0.3,
        nfreq = 30,
        rho = 1000.0,
        c0 = 1480.0,
        aa = 1.0))