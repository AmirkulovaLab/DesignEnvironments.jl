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