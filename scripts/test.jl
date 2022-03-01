using DesignEnvironments
using Plots

## constructing the design
design = Configuration(
    M = 20,
    plane_size = 15.0,
    max_vel = 0.2,
    vel_decay = 0.8,
    min_distance = 0.1)

objective = TSCS(
    k0amax = 2.0,
    k0amin = 0.35,
    nfreq = 100,
    rho = DE.RHO,
    c0 = DE.C0,
    aa = maximum(design.radii))

## locations
x = [
    11.6955    12.9783;
    -12.8262    10.678;
    11.1696    -2.6522;
    -6.01489    1.86352;
    1.74465  -10.0585]

for _ in 1:10
    @time objective(design.pos)
end