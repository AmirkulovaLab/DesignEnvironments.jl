function DE.render(
        env::DesignEnvironment{CoreConfiguration, TSCS}, 
        policy::AbstractPolicy, 
        path::String)

    design = env.design
    tscs = env.objective

    reset!(env)
    freqv = range(tscs.k0amin, tscs.k0amax, length=tscs.nfreq) |> collect

    a = Animation()

    objective_scale = DE.scale(tscs)

    env.objective(design.core)
    core_Q = tscs.Q

    while !is_terminated(env)
        env(policy(env))

        p = img(env, objective_scale)
        plot!(p[2], freqv, core_Q)
        frame(a, p)
    end

    gif(a, path, fps=20)
    closeall()
end