export Render

mutable struct Render <: AbstractHook
    path::String
    a::Animation
    objective_scale::Tuple
end

function Render(path::String)
    return Render(path, Animation(), ())
end

function (render::Render)(::PreEpisodeStage, policy::AbstractPolicy, env::DesignEnvironment)
    reset!(env)
    render.objective_scale = scale(env.objective)
end

function (render::Render)(::PostActStage, policy::AbstractPolicy, env::DesignEnvironment)
    frame(render.a, img(env, render.objective_scale))
end

function (render::Render)(
        ::PostActStage, 
        policy::AbstractPolicy, 
        env::DesignEnvironment{CoreConfiguration, TSCS})

    p = img(env, render.objective_scale)

    k0amin = env.objective.k0amin
    k0amax = env.objective.k0amax
    nfreq = env.objective.nfreq

    freqv = range(k0amin, k0amax, length=nfreq) |> collect
    env.objective(env.design.core)
    plot!(p[2], freqv, env.objective.Q)
    
    frame(render.a, p)
end

function (render::Render)(::PostEpisodeStage, policy::AbstractPolicy, env::DesignEnvironment)
    gif(render.a, render.path, fps=20)
    closeall()
end