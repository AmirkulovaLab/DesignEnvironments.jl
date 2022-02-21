using DesignEnvironments
using LinearAlgebra
using SpecialFunctions: besselj, besselh
using BlockDiagonals: BlockDiagonal
using Plots
import Flux

const RHO = 1000.0
const C0 = 1480.0

mutable struct VectorizedTSCS <: AbstractObjective
    ##
    device::Function

    ## 
    k0amax::Real
    k0amin::Real
    nfreq::Int
    a::Real
    rho::Real
    c0::Real

    ##
    Q::Vector{Real}
    qV::Vector{Real}
end

function VectorizedTSCS(;k0amax::Real, k0amin::Real, nfreq::Int, a::Real, rho::Real, c0::Real)

    if CUDA.has_cuda_gpu()
        device = CuArray
    else
        device = Flux.cpu
    end
    
    return VectorizedTSCS(device, k0amax, k0amin, nfreq, a, rho, c0, [], [])
end

function cartesian_to_polar(x::Vector, y::Vector)
    distances = sqrt.(x .^2 .+ y .^2)
    angles = atan.(y, x)

    return hcat(distances, angles)
end

"""
Calculates the vector of frequencies to calculate tscs at.

# Arguments
-`pa::TSCS`: struct which specifies tscs calculations
"""
function frequency_vector(tscs::VectorizedTSCS)
    freq_max = tscs.k0amax * tscs.c0 / (2 * pi * tscs.a)
    freq_min = tscs.k0amin * tscs.c0 / (2 * pi * tscs.a)

    df = (freq_max - freq_min) / (tscs.nfreq - 1)

    freqv = collect(freq_min:df:freq_max)

    return freqv
end

function t_matrix(nv_i::Matrix, ka_i::Float64)
    J_pv_ka = (besselj.(nv_i .-1, ka_i) - besselj.(nv_i .+ 1, ka_i)) / 2
    H_pv_ka = (besselh.(nv_i .-1, ka_i) - besselh.(nv_i .+ 1, ka_i)) / 2
    return Diagonal(- J_pv_ka ./ H_pv_ka)
end

function plane_wave_amplitude(k0_i::Float64, xM::Vector, nv_i::Matrix, )

end

function (tscs::VectorizedTSCS)(x::Matrix)
    xM = x[:, 1]
    yM = x[:, 2]
    M = length(xM)

    polar_coords = cartesian_to_polar(xM, yM) |> tscs.device
    freqv = frequency_vector(tscs) |> tscs.device

    ## rotating angle by 2 pi if negative
    invalid_angle = polar_coords[:, 2] .< 0
    polar_coords[invalid_angle, 2] .+= 2 * pi

    ## computing scattering at each frequency
    omega = 2 * pi .* freqv
    k0 = omega ./ tscs.c0

    ka = k0 .* tscs.a
    kaa = k0 .* tscs.a ## replace with aa

    nmax = round.(2.5 * ka)

    nv = range.(-nmax, nmax) .|> collect .|> transpose .|> Array

    T0 = t_matrix.(nv, ka)
    T1 = t_matrix.(nv, kaa)

    T = [T0]
    for _ in 1:(M-1)
        push!(T, T1)
    end

    T_diag = BlockDiagonal.(T)

    # exp.(im * k0 * xM[1]) * exp.(im * nv * pi / 2)

    # im * k0 * xM[1]

    k0[1]
end

## locations
x = [
    11.6955    12.9783;
    -12.8262    10.678;
    11.1696    -2.6522;
    -6.01489    1.86352;
    1.74465  -10.0585]

## constructing the design
design = Configuration(
    M = size(x, 1),
    plane_size = 15.0,
    max_vel = 0.2,
    vel_decay = 0.8,
    min_distance = 0.1)

tscs = VectorizedTSCS(
    k0amax = 1.0,
    k0amin = 0.3,
    nfreq = 20,
    a = maximum(design.radii),
    rho = RHO,
    c0 = C0)

tscs(x)