export TSCS

const RHO = 1000.0
const C0 = 1480.0
const aa = 1.0

"""
Calculates the Total Scattering Cross Section of a planar `Configuration` of cylindrical scatterers.
TSCS is calculated independently for each wavenumber (nfreq). Therefore this calculation benifits from
multithreading. Start julia with threads <= nfreq to take advantage of this.

# Example
```
objective = TSCS(k0amax = 1.0, k0amin = 0.3, nfreq = 30)
```

# Parameters
- `k0amax::Float64`
- `k0amin::Float64`
- `nfreq::Int`
- `rho::Float64`
- `c0::Float64`
- `aa::Float64`
- `a::Float64`
- `ha::Float64`
- `Q_RMS::Float64`
- `qV::Vector`
- `Q::Vector`
"""
mutable struct TSCS <: AbstractObjective
    k0amax::Float64
    k0amin::Float64
    nfreq::Int
    rho::Float64
    c0::Float64
    aa::Float64
    a::Float64
    ha::Float64

    ## cache of last calculated objective
    Q_RMS::Float64
    qV::Vector
    Q::Vector
end

function TSCS(; k0amax, k0amin, nfreq, rho, c0, aa)
    a = aa
    ha = aa/10
    Q_RMS = 0.0
    qV = Vector()
    Q = Vector()

    return TSCS(k0amax, k0amin, nfreq, rho, c0, aa, a, ha, Q_RMS, qV, Q)
end

function σ(tscs::TSCS, x::Matrix, xM::Vector, absr::Vector, argr::Vector, freq::Float64)
    ## number of cylinders
    M = size(x, 1)

    ## establishing properties for σ calculation
    omega = 2 * pi * freq
    k0 = omega / tscs.c0
    ka = k0 * tscs.a
    kaa = k0 * tscs.aa
    nmax = round(2.5 * ka)
    N = Int(nmax)
    nv = collect(-nmax:nmax)

    ## T matrix for inner rigid cilinders
    Jpvka = (besselj.(nv .- 1, ka) - besselj.(nv .+ 1, ka)) ./ 2
    Hpvka = (besselh.(nv .- 1, ka) - besselh.(nv .+ 1, ka)) ./ 2
    T_0 = Diagonal(- Jpvka ./ Hpvka)

    ## T matrix for outer cloaking cilinders
    Jpvkaa = (besselj.(nv .- 1, kaa) - besselj.(nv .+ 1, kaa)) ./ 2
    Hpvkaa = (besselh.(nv .- 1, kaa) - besselh.(nv .+ 1, kaa)) ./ 2
    T_1 = Diagonal(- Jpvkaa ./ Hpvkaa)

    ## creating 3d matrix
    T = [T_0]
    for _ in 2 : M
        append!(T, [T_1])
    end

    ## creating big block diagonal matrix
    Tdiag = BlockDiagonal(T)

    N_size = 2 * N + 1

    Ainv = zeros(ComplexF64, N_size, M)
    AM = zeros(ComplexF64, N_size, M)

    ## incident plane wave amplitude
    Ainv[:, 1] = exp(1im * k0 * xM[1]) * exp.(nv .* (1im * pi / 2))
    AM[:, 1] = T_0 * Ainv[:, 1]

    ## incident plane wave amplitude
    Ainv[:, 1] = exp(1im * k0 * xM[1]) * exp.(nv .* (1im * pi / 2))
    AM[:, 1] = T_0 * Ainv[:, 1]

    for j = 2 : M
        Ainv[:, j] = exp(1im * k0 * xM[j]) * exp.(nv .* (1im * pi / 2))
        AM[:, j] = T_1 * Ainv[:, j]
    end

    verAv = reshape(AM, M * N_size, 1)

    Xbig = Matrix{ComplexF64}(I, M * N_size, M * N_size)
    fill!(Xbig, 0 + 0im)

    ## defining tensors to hold gradients
    absrjm = zeros(M, M)
    argrjm = zeros(M, M)
    P_jm = zeros(ComplexF64, N_size, N_size)
    gradP_jm = zeros(ComplexF64, N_size, N_size, M, M, 2)
    gradP_mj = zeros(ComplexF64, N_size, N_size, M, M, 2)
    gradXbigM = zeros(ComplexF64, N_size, N_size, M, M, 2)
    gradXbigM1 = zeros(ComplexF64, N_size, N_size, M, M, 2)
    gradXbig = zeros(ComplexF64, M * N_size, M * N_size, M, 2)
    Dp_jm = zeros(ComplexF64, N_size, N_size)
    Dh_jm = zeros(ComplexF64, N_size, N_size)
    Dp_mj = zeros(ComplexF64, N_size, N_size)
    Dh_mj = zeros(ComplexF64, N_size, N_size)

    if M == 1
        ## edge case for single cylinder
        Xbig = Matrix{ComplexF64}(I, N_size, N_size)
        gradXbig = zeros(ComplexF64, N_size, N_size, M, 2)
    else
        ## case for cylinders > 1
        for j = 1 : M
            for m = 1 : M
                j_idx = ((j - 1) * N_size + 1 : j * N_size)
                m_idx = ((m - 1) * N_size + 1 : m * N_size)

                if m == j
                    ## setting diagonal of big matrix to identity
                    Xbig[j_idx, m_idx] = I(N_size)
                else
                    rjm = x[j, :] - x[m, :]
                    xjm = rjm[1]
                    yjm = rjm[2]

                    ## distance between two cylinders
                    absrjm[j, m] = norm(rjm)
                    argrjm[j, m] = atan(yjm, xjm)

                    ## shifting by two periods if negative
                    if argrjm[j, m] < 0
                        argrjm[j, m] += 2 * pi
                    end

                    for i = 1 : length(nv)
                        n = nv[i]
                        Hvrjm = besselh.(n .- nv, k0 * absrjm[j, m])
                        Hpvrjm = (besselh.(n .- nv .- 1, k0 * absrjm[j, m]) - besselh.(n .- nv .+ 1, k0 * absrjm[j, m])) ./ 2
                        exprjm = exp.(1im .* (n .- nv) .* argrjm[j, m])

                        Dp_jm[:, i] = Hpvrjm .* exprjm
                        P_jm[:, i] = Hvrjm .* exprjm
                        Dh_jm[:, i] = 1im .* (n .- nv) .* Hvrjm .* exprjm

                        exprmj = exp.(1im .* (n .- nv) .* (argrjm[j, m] .+ pi))
                        Dp_mj[:, i] = Hpvrjm .* exprmj
                        Dh_mj[:, i] = 1im .* (n .- nv) .* Hvrjm .* exprmj
                    end

                    Xbig[j_idx, m_idx] = -Tdiag[j_idx, j_idx] * P_jm

                    gradP_jm[:, :, m, j, 1] = (k0 * xjm * Dp_jm / absrjm[j, m]) - (Dh_jm * yjm / absrjm[j, m]^2)
                    gradP_jm[:, :, m, j, 2] = (k0 * yjm * Dp_jm / absrjm[j, m]) + (Dh_jm * xjm / absrjm[j, m]^2)

                    gradXbigM[:, :, m, j, 1] = -Tdiag[j_idx, j_idx] * gradP_jm[:, :, m, j, 1]
                    gradXbigM[:, :, m, j, 2] = -Tdiag[j_idx, j_idx] * gradP_jm[:, :, m, j, 2]

                    gradXbig[j_idx, m_idx, j, 1] = gradXbigM[:, :, m, j, 1]
                    gradXbig[j_idx, m_idx, j, 2] = gradXbigM[:, :, m, j, 2]

                    gradP_mj[:, :, m, j, 1] = (k0 * xjm * Dp_mj / absrjm[j, m]) - (Dh_mj * yjm / absrjm[j, m]^2)
                    gradP_mj[:, :, m, j, 2] = (k0 * yjm * Dp_mj / absrjm[j, m]) + (Dh_mj * xjm / absrjm[j, m]^2)

                    gradXbigM1[:, :, m, j, 1] = -Tdiag[j_idx, j_idx] * gradP_mj[:, :, m, j, 1]
                    gradXbigM1[:, :, m, j, 2] = -Tdiag[j_idx, j_idx] * gradP_mj[:, :, m, j, 2]

                    gradXbig[m_idx, j_idx, j, 1] = gradXbigM1[:, :, m, j, 1]
                    gradXbig[m_idx, j_idx, j, 2] = gradXbigM1[:, :, m, j, 2]
                end
            end
        end
    end

    ## solving linear system
    Bv = Xbig \ verAv
    Bcoef = reshape(Bv, N_size, M)

    C_nm = zeros(ComplexF64, N_size, M)

    for i = 1 : length(nv)
        n = nv[i]
        for m = 1 : M
            C_nm[i, m] = (-1im)^n * exp(-1im * k0 * absr[m] * cos(argr[m]))
        end
    end

    Cv = reshape(C_nm, M * N_size, 1)

    ## filling gradient matrix with complex zeros
    gradcM = zeros(ComplexF64, N_size, M, M, 2)
    GradAinv = zeros(ComplexF64, N_size, M, M, 2)
    gradaM = zeros(ComplexF64, N_size, M, M, 2)
    gradAv = zeros(ComplexF64, N_size * M, M, 2)
    gradCv = zeros(ComplexF64, N_size * M, M, 2)

    for j = 1 : M
        for m = 1 : M
            if m == j
                gradcM[:, m, j, 1] = -1im * k0 * C_nm[:, m]
                gradcM[:, m, j, 2] .= 0

                GradAinv[:, m, j, 1] = 1im * k0 * Ainv[:, m]
                GradAinv[:, m, j, 2] .= 0
            end
        end
    end

    for i = 1 : 2
        for j = 1 : M
            j_idx = ((j - 1) * N_size + 1 : j * N_size)

            gradaM[:, :, j, i] = Tdiag[j_idx, j_idx] * GradAinv[:, :, j, i]
            gradAv[:, j, i] = reshape(gradaM[:, :, j, i], M * N_size, 1)
            gradCv[:, j, i] = reshape(gradcM[:, :, j, i], M * N_size, 1)
        end
    end

    ## computing tscs
    Q = (-4 / k0 * real(transpose(Cv) * Bv))[1]

    s_j = zeros(M, 2)
    q_j = zeros(M, 2)

    ## computing gradient of tscs
    for i = 1 : 2
        for j = 1 : M
            left = transpose(gradCv[:, j, i]) * Bv
            right = transpose(Cv) * (Xbig \ (gradXbig[:, :, j, i] * Bv - gradAv[:, j, i]))

            s_j[j, i] = - 4 / k0 * real(left - right)[1]
            q_j[j, i] = Q * s_j[j, i]
        end
    end

    return (Q, q_j)
end

"""
Defines the function call for `TSCS` on a `Matrix` of scatterer coordinates.

# Example
```
objective = TSCS(k0amax = 1.0, k0amin = 0.3, nfreq = 30)
objective(x)
```
"""
function (tscs::TSCS)(x::Matrix)
    M = size(x, 1)

    ## reshaping and transpose
    # x = Matrix(reshape(x, 2, M)')
    xM = x[:,1]
    yM = x[:,2]

    absr = sqrt.((xM .^ 2) + (yM .^ 2))
    argr = atan.(yM, xM)

    ## shifting negative values by 2pi
    argr[argr .< 0] .+= (2 * pi)

    # getting frequency range
    freqmax = (tscs.k0amax * tscs.c0) / (2 * pi * tscs.a)
    freqmin = (tscs.k0amin * tscs.c0) / (2 * pi * tscs.a)
    df = (freqmax - freqmin) / (tscs.nfreq - 1)

    freqv = collect(range(freqmin, freqmax, step=df))

    Q = zeros(tscs.nfreq)
    q_j = zeros(tscs.nfreq, M, 2)

    Threads.@threads for Ifreq = 1 : length(freqv)
        freq = freqv[Ifreq][1]
        Q[Ifreq], q_j[Ifreq, :, :] = σ(tscs, x, xM, absr, argr, freq)
    end

    ## computing final values
    Q_RMS = sqrt((1 / tscs.nfreq) * (sum(Q .^ 2)))
    q_RMS_j = 1 / (tscs.nfreq * Q_RMS) .* sum(q_j, dims=1)
    q_RMS_j = dropdims(q_RMS_j, dims=1)
    qV = reshape(transpose(q_RMS_j), 2 * M, 1)
    qV = Vector{Float64}(vec(qV))

    ## setting cache
    tscs.Q_RMS = Q_RMS
    tscs.qV = qV
    tscs.Q = Q
end

"""
Defines function call for `TSCS` on a `Configuration`

# Example
```
objective(config)
```
"""
function (tscs::TSCS)(config::Configuration)
    return tscs(config.pos)
end

"""
Defines function call for `TSCS` on a `CoreConfiguration`
"""
function (tscs::TSCS)(design::CoreConfiguration)
    config = merge_configs(design.core, design.config)
    return tscs(config)
end

function scale(tscs::TSCS)
    return (0.0, maximum(tscs.Q))
end

function metric(tscs::TSCS)
    return tscs.Q_RMS
end

function Plots.plot(tscs::TSCS, objective_scale::Tuple)
    freqv = range(tscs.k0amin, tscs.k0amax, length=tscs.nfreq) |> collect

    return plot(
        freqv, tscs.Q,
        xlabel="ka", ylabel="TSCS",
        xlim = (freqv[1], freqv[end]),
        ylim = objective_scale, legend=false)
end

function Plots.plot(tscs::TSCS)
    return plot(tscs, scale(tscs))
end

