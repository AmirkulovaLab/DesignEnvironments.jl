using DesignEnvironments
using Plots
using LinearAlgebra

const RHO = 1000.0
const C0 = 1480.0

struct PressureAmplitude <: AbstractObjective
    k0amax::Real
    k0amin::Real
    nfreq::Int
    R2::Real
    a::Real
    rho::Real
    c0::Real
end

function cartesian_to_polar(x::Vector, y::Vector)
    distances = sqrt.(x .^2 .+ y .^2)
    angles = atan.(y, x)

    return hcat(distances, angles)
end

"""
Calculates relevent geometrical relationships between cylinders, far field point, and focal point.

# Arguments
-`xM::Vector`: x coords of cylinders
-`yM::Vector`: y coords of cylinders
-`radii::Vector`: radii of all cylinders
-`focal_x::Real`: x coord of focal point
-`focal_y::Real`: y coord of focal point
"""
function geometry(
        xM::Vector, yM::Vector, a::Real,
        focal_x::Real, focal_y::Real)

    ## obtaining distances and angles of all cylinders from origin
    distances = sqrt.(xM .^2 .+ yM .^2)
    angles = atan.(yM, xM)
    cylinder_polar = hcat(distances, angles)

    ## obtaining distances and angles of cylinders relative to focal point
    cylinder_focal_distance = sqrt.((focal_x .- xM) .^2 + (focal_y .- yM) .^2)
    cylinder_focal_angle = atan.((focal_y .- yM), (focal_x .- xM))
    cylinder_focal_polar = hcat(cylinder_focal_distance, cylinder_focal_angle)

    ## computing far field point based on largest cylinder radii
    far_field_point = [-100 * a, 0.0]
    far_field_distance = norm(far_field_point)

    return cylinder_polar, cylinder_focal_polar, far_field_distance
end

"""
Calculates the vector of frequencies to calculate pressure amplitude at.

# Arguments
-`pa::PressureAmplitude`: struct which specifies pressure amplitude calculations
"""
function frequency(pa::PressureAmplitude)
    freq_max = pa.k0amax * pa.c0 / (2 * pi * pa.a)
    freq_min = pa.k0amin * pa.c0 / (2 * pi * pa.a)

    df = (freq_max - freq_min) / (pa.nfreq - 1)

    freqv = collect(freq_min:df:freq_max)

    return freqv
end

function (pa::PressureAmplitude)(x::Matrix,  xf::Vector)
    xM = x[:, 1] ## x coords of cylinders
    yM = x[:, 2] ## y coords of cylinders
    focal_x, focal_y = xf ## x and y coords of focal point
    M = length(xM) ## obtaining number of cylinders

    ## calculating necissary geometric relationships
    cylinder_polar, cylinder_focal_polar, far_field_distance = 
        geometry(xM, yM, pa.a, focal_x, focal_y)

    ## obtaining frequency vector for PressureAmplitude
    freqv = frequency(pa)

    ## initializing zero arrays to hold data
    kav = zeros(1, size(freqv)...)
    Q = zeros(size(freqv)... , 1)
    q_j = zeros(pa.nfreq, M, 2)
    q_j2 = deepcopy(q_j)

    ## rotating frequency by one period
    omega = 2 * pi .* freqv
    k0 = omega ./ pa.c0

end

## constructing the design
design = Configuration(
    M = 5,
    plane_size = 15.0,
    max_vel = 0.2,
    vel_decay = 0.8,
    min_distance = 0.1)

pa = PressureAmplitude(
    0.45, 
    0.35, 
    11, 
    10.0, 
    maximum(design.radii),
    RHO, 
    C0)

## focal point
xf = [12.0, 0.0]

## locations
x = [
    11.6955    12.9783;
    -12.8262    10.678;
    11.1696    -2.6522;
    -6.01489    1.86352;
    1.74465  -10.0585]

pa(x, xf)