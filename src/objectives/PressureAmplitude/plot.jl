function Plots.plot(pa::PressureAmplitude, objective_scale::Tuple)
    freqv = collect(range(pa.k0amin, pa.k0amax, length=pa.nfreq))

    return plot(
        freqv, pa.Q,
        xlabel = "ka", ylabel = "Pressure Amplitude",
        xlim = (freqv[1], freqv[end]),
        ylim = objective_scale, 
        legend=false)
end

function Plots.plot(pa::PressureAmplitude)
    return plot(pa, scale(pa))
end