using DesignEnvironments

env = CylinderEnv(
    M=10, 
    plane=Square(10.0),
    nfreq=20,
    k0amax=1.0,
    max_vel=0.1,
    vel_decay=1.0)

