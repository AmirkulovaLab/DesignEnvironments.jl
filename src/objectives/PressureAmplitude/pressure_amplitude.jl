export PressureAmplitude

mutable struct PressureAmplitude <: AbstractObjective
    ## focal point
    xf::AbstractArray

    use_cuda::Bool
    k0amax::Real
    k0amin::Real
    nfreq::Int
    R2::Real
    a::Real
    rho::Real
    c0::Real

    ## last calculated Q
    Q::Vector
end

function PressureAmplitude(;
        xf, use_cuda, k0amax, k0amin, nfreq, R2, a, rho, c0)
    Q = zeros(nfreq)

    return PressureAmplitude(xf, use_cuda,  k0amax,  k0amin,  nfreq, R2, a, rho, c0, Q)
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
Construct a matrix with the inverse T matrix repeated along the diagonal M times.
"""
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

"""
Obtain square matricies containing distances and angles between each cylinder. Zero along the diagonal.
"""
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

function fill_Xbig!(Xbig::Matrix, nv::Vector, nmax::Int, k0::Real, absrjm::Matrix, argrjm::Matrix)

    M = size(absrjm, 1)
    Dp_jm = Vector{Vector{ComplexF32}}()
    Dp_mj = Vector{Vector{ComplexF32}}()
    Dh_jm = Vector{Vector{ComplexF32}}()
    Dh_mj = Vector{Vector{ComplexF32}}()
    P_jm = Vector{Vector{ComplexF32}}()

    for j = 1:M
        for m = 1:M
            if j != m
                P_jm = Vector{Vector{ComplexF32}}()
                scaled_distance = k0 * absrjm[j, m]
                for n in nv
                    shifted_nv = n .- nv
                    Hvrjm = besselh.(shifted_nv, scaled_distance)
                    Hpvrjm = (besselh.(shifted_nv .- 1, scaled_distance) .- besselh.(shifted_nv .+ 1, scaled_distance)) / 2
                    exprjm = exp.(im * shifted_nv * argrjm[j, m])

                    push!(Dp_jm, Hpvrjm .* exprjm)

                    Hv_exp = Hvrjm .* exprjm

                    push!(P_jm, Hv_exp)

                    scaled_Hv_exp = im * shifted_nv .* Hv_exp
                    push!(Dh_jm, scaled_Hv_exp)

                    exprmj = exp.(im * shifted_nv * (argrjm[j, m] + pi))

                    push!(Dp_mj, Hpvrjm .* exprmj)
                    push!(Dh_mj, im * shifted_nv .* Hvrjm .* exprmj)
                end
                P_jm = hcat(P_jm...)

                N = nmax
                j_start = (j - 1) * (2 * N + 1) + 1
                j_end = j * (2 * N + 1)

                m_start = (m - 1) * (2 * N + 1) + 1
                m_end = m * (2 * N + 1)
                Xbig[j_start:j_end, m_start:m_end] = -P_jm
            end
        end
    end

    solution_matricies = (Dp_jm, Dp_mj, Dh_jm, Dh_mj, P_jm)
    solution_matricies = (hcat(solution...) for solution in solution_matricies)
    return solution_matricies
end


function focal_point_calculation(x::Matrix, xf::Vector, nv::Vector, k0::Real)
    M = size(x, 1)
    N = Int(maximum(nv))

    Vh_fm= zeros(ComplexF32, 2*N+1, M)

    for j = 1:M
        for m = 1:M
            rfm = xf - x[m, :]
            xfm, yfm = rfm
            absrfm = norm(rfm)
            argrfm = atan(yfm, xfm)

            exprfm = exp.(im * nv * argrfm)
            
            scaled_angle = k0 * absrfm
            Hvrfm = besselh.(nv, scaled_angle)
            Jvrfm = besselj.(nv, scaled_angle)

            Hpvrfm = (besselh.(nv .- 1, scaled_angle) .- besselh.(nv .+ 1, scaled_angle)) / 2
            Jpvrfm = (besselj.(nv .- 1, scaled_angle) .- besselj.(nv .+ 1, scaled_angle)) / 2

            Dp_fm = Hpvrfm .* exprfm
            Dh_fm = im * nv .* Hvrfm .* exprfm
            Vh_fm[:, m] = Hvrfm .* exprfm;
        end
    end

    return Vh_fm
end

"""
main function which is called to compute pressure at focal point xf produced by configuration of scatterers x
"""
function (pa::PressureAmplitude)(x::Matrix)
    xf = pa.xf

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
    nmax = Vector{Float64}(round.(Int, 2.5 * ka, RoundNearestTiesUp))
    N = Int.(nmax)
    nv = collect.(range.(-nmax, nmax))

    ## incident plane wave
    ## normal incident wave direction
    psi_0deg = 0
    psi_0 = psi_0deg * pi / 180 
    p_i = exp.(im * k0 * (xf[1] * cos(psi_0) + xf[2] * sin(psi_0)))

    ## comuting pressure vector for all ka
    pv_ka = pressure_vector.(nv, ka)
    J_pv_ka = getindex.(pv_ka, 1)
    H_pv_ka = getindex.(pv_ka, 2)

    ## computing t matrix
    # T0 = t_matrix.(J_pv_ka, H_pv_ka)
    # T1 = deepcopy(T0)
    invT_0 = t_matrix.(H_pv_ka, J_pv_ka)

    # D0 = sqrt.(im * 0.5 * pi * k0 * far_field_distance) .* exp.(-im * k0 * far_field_distance)

    Ainv = []
    for j = 1 : M
        ## create a matrix of (nv, nfreq) for each cylinder j
        Ainv_i = exp.(im * k0 * xM[j]) .* im_vector_power.(nv)
        push!(Ainv, Ainv_i)
    end

    # ## get dimention of flattened matrix
    vec_dim = (M * (2 * N .+ 1))

    # Ainv = transpose.(hcat.(Ainv...))
    Ainv = hcat.(Ainv...)

    ## obtain a vector of flattened matricies (vectors)
    verAv = Vector.(reshape.(Ainv, vec_dim))
    Xbig = build_Xbig.(invT_0, M)

    ## computing distance and angle matricies
    rjm = build_rjm(x)
    absrjm = rjm[1, :, :]
    argrjm = rjm[2, :, :]

    fill_Xbig!.(
        Xbig,
        nv,
        N,
        k0,
        [absrjm for _ in 1:pa.nfreq],
        [argrjm for _ in 1:pa.nfreq])

    if pa.use_cuda
        verAv = cu.(verAv)
        Xbig = cu.(Xbig)
    else
        verAv = Vector{ComplexF32}.(verAv)
        Xbig = Matrix{ComplexF32}.(Xbig)
    end

    bV = Xbig .\ verAv

    Vh_fm = focal_point_calculation.(
        [x for _ in 1:pa.nfreq],
        [xf for _ in 1:pa.nfreq],
        nv, 
        k0)

    Vv_fm = transpose.(reshape.(Vh_fm, M * (2 * N .+ 1)))

    if pa.use_cuda
        Vv_fm = cu.(Vv_fm)
        p_i = cu.(p_i)
    end

    pf_fm1 = p_i .+ Vv_fm .* bV
    pa.Q = abs.(pf_fm1)
end

function (pa::PressureAmplitude)(x::Matrix, xf::Matrix)
    pa.(x, xf)
end

function scale(pa::PressureAmplitude)
    return (0.0, maximum(pa.Q))
end

function metric(pa::PressureAmplitude)
    return sqrt(mean(pa.Q .^ 2))
end