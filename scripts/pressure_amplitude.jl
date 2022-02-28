using DesignEnvironments
using Plots
using LinearAlgebra
using SpecialFunctions: besselh, besselj
using BlockDiagonals: BlockDiagonal
using CUDA

const RHO = 1000.0
const C0 = 1480.0

Base.@kwdef struct PressureAmplitude <: AbstractObjective
    use_cuda::Bool
    k0amax::Real
    k0amin::Real
    nfreq::Int
    R2::Real
    a::Real
    rho::Real
    c0::Real
end

"""
Convert vectors of x and y coordinates into polar coordinate matrix with columns (radius, angle)
"""
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

"""
Computes J and H component of the pressure vector over a range of nv.
"""
function pressure_vector(nv::Vector, ka::Real)
    J_pv_ka = (besselj.(nv' .- 1, ka) - besselj.(nv' .+ 1, ka)) / 2
    H_pv_ka = (besselh.(nv' .-1, ka) - besselh.(nv' .+ 1, ka)) / 2
    return vec(J_pv_ka), vec(H_pv_ka)
end

"""
Take in J_pv_ka and H_pv_ka vectors for single frequency. Solve linear system and place it as diagonal
in a zero matrix.
"""
function t_matrix(J_pv_ka::Vector, H_pv_ka::Vector)
    return - (J_pv_ka ./ H_pv_ka) |> diagm
end

"""
raise the imaginary unit to the power of a vector.
"""
function im_vector_power(nv_i::Vector)
    return im .^ nv_i
end

"""
main function which is called to compute pressure at focal point xf produced by configuration of scatterers x
"""
function (pa::PressureAmplitude)(x::Matrix,  xf::Vector)
    xM = x[:, 1] ## x coords of cylinders
    yM = x[:, 2] ## y coords of cylinders
    focal_x, focal_y = xf ## x and y coords of focal point
    M = length(xM) ## obtaining number of cylinders

    ## calculating necissary geometric relationships
    cylinder_polar, cylinder_focal_polar, far_field_distance = geometry(xM, yM, pa.a, focal_x, focal_y)

    ## obtaining frequency vector for PressureAmplitude
    freqv = frequency(pa)

    ## rotating frequency by one period
    omega = 2 * pi .* freqv
    k0 = omega ./ pa.c0

    ka = k0 .* pa.a
    kaa = k0 .* pa.a ## replace with aa
    nmax = round.(2.5 * ka)
    nv = range.(-nmax, nmax) .|> collect

    ## comuting pressure vector for all ka
    pv_ka = pressure_vector.(nv, ka)
    J_pv_ka = getindex.(pv_ka, 1)
    H_pv_ka = getindex.(pv_ka, 2)

    ## computing t matrix
    T0 = t_matrix.(J_pv_ka, H_pv_ka)
    T1 = deepcopy(T0)
    invT_0 = t_matrix.(H_pv_ka, J_pv_ka)
    # D0 = sqrt.(im * 0.5 * pi * k0 * far_field_distance) .* exp.(-im * k0 * far_field_distance)

    Ainv = []
    for j = 1 : M
        ## create a matrix of (nv, nfreq) for each cylinder j
        Ainv_i = exp.(im * k0 * xM[j]) .* im_vector_power.(nv)
        ## concatenate vector of vectors into matrix
        Ainv_i = hcat(Ainv_i...)
        ## store matrix into vector
        push!(Ainv, Ainv_i)
    end

    ## stack matricies along new dimention
    Ainv = cat(Ainv..., dims = 3) |> cu

    ## extract matricies along the frequency dimention
    Ainv = [Ainv[:, i, :] for i in 1:pa.nfreq]
    ## get dimention of flattened matrix
    vec_dim = Int.(M * (2 * nmax .+ 1))
    ## obtain a vector of flattened matricies (vectors)
    verAv = vec.(reshape.(Ainv, vec_dim, 1))

    hcat(verAv...)
end

## constructing the design
design = Configuration(
    M = 5,
    plane_size = 15.0,
    max_vel = 0.2,
    vel_decay = 0.8,
    min_distance = 0.1)

pa = PressureAmplitude(
    use_cuda = false,
    k0amax = 0.45, 
    k0amin = 0.35, 
    nfreq = 11, 
    R2 = 10.0, 
    a = maximum(design.radii),
    rho = RHO, 
    c0 = C0)

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