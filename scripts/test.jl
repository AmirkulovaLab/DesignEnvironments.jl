using DesignEnvironments

design_params = Dict(
    :M => 20,
    :plane_size => 50.0,
    :max_vel => 0.2,
    :vel_decay => 0.8,
    :min_distance => 0.1)

## constructing the design
design = Configuration(; design_params...)

objective_params = Dict(
    :k0amax => 2.0,
    :k0amin => 0.35,
    :nfreq => 30,
    :R2 => design_params[:plane_size],
    :a => maximum(design.radii),
    :rho => DE.RHO,
    :c0 => DE.C0)

pa = PressureAmplitude(
    use_cuda = false;
    objective_params...)

## focal point
xf = [12.0, 0.0]

## locations
x = [
    11.6955    12.9783;
    -12.8262    10.678;
    11.1696    -2.6522;
    -6.01489    1.86352;
    1.74465  -10.0585]

for _ in 1:10
    DE.reset_design!(design)
    @time pa(design.pos, xf)
end