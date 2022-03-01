function Plots.plot(design::AbstractDesign, objective::AbstractObjective, objective_scale::Tuple)
    return plot(
        plot(design), 
        plot(objective, objective_scale), 
        layout=@layout([a{0.6w} b]))
end

function Plots.plot(design::AbstractDesign, objective::AbstractObjective)
    return plot(design, objective, scale(objective))
end

function Plots.plot(env::DesignEnvironment, objective_scale::Tuple)
    return plot(env.design, env.objective, objective_scale)
end

function Plots.plot(env::DesignEnvironment)
    return plot(env, scale(env.objective))
end