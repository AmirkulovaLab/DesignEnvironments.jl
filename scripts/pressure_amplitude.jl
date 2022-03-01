using DesignEnvironments
using Plots
using LinearAlgebra
using SpecialFunctions: besselh, besselj
using BlockDiagonals: BlockDiagonal
using SparseArrays: blockdiag, sparse
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

function build_Xbig(invT::Matrix, M::Int)
    vec_dim = size(invT, 1)
    Xbig = zeros(ComplexF64, vec_dim * M, vec_dim * M)

    for i = 1:M
        idx = vec_dim * (i - 1) + 1
        idx_end = vec_dim * i
        Xbig[idx:idx_end, idx:idx_end] = invT
    end

    return Xbig
end

function build_rjm(x::Matrix)
    M = size(x, 1)
    rjm = zeros(2, M, M)

    for j = 1:M
        for m = 1:M
            if j != m
                distance = x[j, :] - x[m, :]
                xjm, yjm = distance
                rjm[1, j, m] = norm(distance)
                rjm[2, j, m] = atan(yjm, xjm)
            end
        end
    end

    return rjm
end

function solve_bessel(nv::Vector, k0::Real, absrjm::Matrix)

    M = size(absrjm, 1)
    Dp_jm = Vector{Vector{ComplexF32}}()
    Dp_mj = Vector{Vector{ComplexF32}}()
    Dh_jm = Vector{Vector{ComplexF32}}()
    Dh_mj = Vector{Vector{ComplexF32}}()
    P_jm = Vector{Vector{ComplexF32}}()

    for j = 1:M
        for m = 1:M
            if j != m
                scaled_distance = k0 * absrjm[j, m]
                for n in nv
                    shifted_nv = n .- nv
                    Hvrjm = besselh.(shifted_nv, scaled_distance)
                    Hpvrjm = besselh.(shifted_nv .- 1, scaled_distance) .- besselh.(shifted_nv .+ 1, scaled_distance) / 2

                    exprjm = exp.(im * shifted_nv * scaled_distance)
                    push!(Dp_jm, Hpvrjm .* exprjm)

                    Hv_exp = Hvrjm .* exprjm
                    push!(P_jm, Hv_exp)

                    scaled_Hv_exp = im * shifted_nv .* Hv_exp
                    push!(Dh_jm, scaled_Hv_exp)

                    exprmj = exp.(im * shifted_nv * (absrjm[j, m] + pi))
                    push!(Dp_mj, Hpvrjm .* exprmj)
                    push!(Dh_mj, im * shifted_nv .* Hvrjm .* exprmj)
                end
            end
        end
    end

    solution_matricies = (Dp_jm, Dp_mj, Dh_jm, Dh_mj, P_jm)
    solution_matricies = (hcat(solution...) for solution in solution_matricies)
    return solution_matricies
end

"""
main function which is called to compute pressure at focal point xf produced by configuration of scatterers x
"""
function (pa::PressureAmplitude)(x::Matrix, xf::Vector)
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
    nv = collect.(range.(-nmax, nmax))

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
        push!(Ainv, Ainv_i)
    end

    # ## get dimention of flattened matrix
    vec_dim = (M * (2 * Int.(nmax) .+ 1))

    Ainv = transpose.(hcat.(Ainv...))
    ## obtain a vector of flattened matricies (vectors)
    verAv = Vector.(reshape.(Ainv, vec_dim))
    Xbig = build_Xbig.(invT_0, M)

    ## computing distance and angle matricies
    rjm = build_rjm(x)
    absrjm = rjm[1, :, :]
    argrjm = rjm[2, :, :]

    # scaled_distances = [k0[i] * absrjm for i = 1:length(k0)]

    solution = solve_bessel.(
        nv,
        k0,
        [absrjm for _ in 1:pa.nfreq])

    # solve_bessel.(nv, repeat(absrjm, pa.nfreq))

    # if pa.use_cuda
    #     verAv = cu.(verAv)
    #     Xbig = cu.(Xbig)
    # else
    #     verAv = Vector{ComplexF32}.(verAv)
    #     Xbig = Matrix{ComplexF32}.(Xbig)
    # end

    # bV = Xbig .\ verAv
    # return bV
end

design_params = Dict(
    :M => 5,
    :plane_size => 50.0,
    :max_vel => 0.2,
    :vel_decay => 0.8,
    :min_distance => 0.1)

## constructing the design
design = Configuration(; design_params...)

objective_params = Dict(
    :k0amax => 0.45,
    :k0amin => 0.35,
    :nfreq => 11,
    :R2 => design_params[:plane_size],
    :a => maximum(design.radii),
    :rho => RHO,
    :c0 => C0)

pa_cpu = PressureAmplitude(
    use_cuda = false;
    objective_params...)

pa_cuda = PressureAmplitude(
    use_cuda = true;
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

test_iterations = 5

display("CPU: ")
for i = 1:test_iterations
    @time pa_cpu(design.pos, xf)
end

display("CUDA: ")
for i = 1:test_iterations
    @time pa_cuda(design.pos, xf)
end