export PressureAmplitude

struct PressureAmplitude <: AbstractObjective
    M::Int
    k0amax::Real
    k0amin::Real
    nfreq::Int
    R2::Real
    a::Real
    aa::Real
end

function (pa::PressureAmplitude)(x::Matrix)
    return x
end

pa = PressureAmplitude(4, 1.0, 0.3, 11, 10.0, 1.0, 1.0)

