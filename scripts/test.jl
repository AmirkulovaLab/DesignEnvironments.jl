using DesignEnvironments
using Plots: plot, savefig

design_params = Dict(
    :M => 10,
    :plane_size => 15.0,
    :max_vel => 0.2,
    :vel_decay => 0.8,
    :min_distance => 0.1)

## constructing the design
cylinders = Configuration(; design_params...)

objective_params = Dict(
    :use_cuda => false,
    :k0amax => 1.0,
    :k0amin => 0.35,
    :nfreq => 50,
    :R2 => design_params[:plane_size],
    :a => maximum(cylinders.radii),
    :rho => DE.RHO,
    :c0 => DE.C0)

pressure_amplitude = PressureAmplitude(
    xf = [12.0, 0.0];
    objective_params...)

env = DesignEnvironment(
    design = cylinders,
    objective = pressure_amplitude,
    is_continuous = false,
    episode_length = 100,
    penalty_weight = 0.1,
    state_type = SequenceVectorState)

savefig(plot(env), "env")

# ## locations
# x = [
#     11.6955    12.9783;
#     -12.8262    10.678;
#     11.1696    -2.6522;
#     -6.01489    1.86352;
#     1.74465  -10.0585]

# for _ in 1:10
#     DE.reset_design!(design)
#     @time pa(design.pos, xf)
# end