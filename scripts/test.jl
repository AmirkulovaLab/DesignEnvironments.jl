using DesignEnvironments
using ReinforcementLearning
using Plots
using IntervalSets

## define objectives

cylinders = Configuration(8, 10.0, 0.2, 0.7, 0.1)
pressure_amplitude = PressureAmplitude(
    xf = [12.0, 0.0],
    k0amax = 0.45,
    k0amin = 0.35,
    nfreq = 10,
    a = maximum(cylinders.radii),
    R2 = cylinders.plane.size,
    rho = DE.RHO,
    c0 = DE.C0,
    use_cuda = false)

env = DesignEnvironment(
    design = cylinders,
    # objective = TSCS(1.0, 0.35, 10, DE.RHO, DE.C0)
    objective = pressure_amplitude,
    episode_length = 100,
    penalty_weight = 0.1,
    is_continuous = true,
    state_type = SequenceState)

size(state_space(env))

# env = DesignEnvironment(
#     design = cylinders,
#     objective = pressure_amplitude,
#     episode_length = 100,
#     penalty_weight = 0.1,
#     is_continuous = true,
#     state_type = VectorState)

# pressure_amplitude(cylinders)

# savefig(plot(env), "pa")

# function focal_ring(x::Real, y::Real, r::Real, n::Int)
#     theta = (2 * pi) / n

#     ring = Vector()

#     for i = 1 : n

#         x_r_i = x + r*cos(theta*i)
#         y_r_i = y + r*sin(theta*i)

#         xy = [x_r_i y_r_i]
#         push!(ring, xy)
#     end

#     return vcat(ring...)
# end

# ring1 = focal_ring(12.0, 0.0, 1.0, 10)
# ring2 = focal_ring(12.0, 0.0, 0.5, 10)

# p = scatter(ring1[:, 1], ring1[:, 2], aspect_ratio = :equal)
# p2 = scatter!(p, ring2[:, 1], ring2[:, 2])
# savefig(p2, "ring")



