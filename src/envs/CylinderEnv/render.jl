export img, render

function cylinder(x, y, r=1; n=30)
    θ = 0:360÷n:360
    Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
end

function img(env::CylinderEnv)

    coords = get_coords(env)
    x, y = coords[:, 1], coords[:, 2]
    cylinders = cylinder.(x, y)

    p = plot(
        cylinders;
        aspect_ratio=:equal,
        legend=false,
        color=:black,
        xlim=(-env.params.grid_size, env.params.grid_size),
        ylim=(-env.params.grid_size, env.params.grid_size)
        )

    return p
end

function render(
        env::CylinderEnv, policy::AbstractPolicy;
        path::String="anim.mp4", fps::Int=30
        )

    k0amin = env.params.k0amin
    k0amax = env.params.k0amax
    nfreq = env.params.nfreq

    dif = (k0amax - k0amin) / nfreq
    freqv = range(k0amin, k0amax, length=nfreq) |> collect

    reset!(env)

    a = Animation()

    while !is_terminated(env)

        p = plot(
            img(env),
            plot(freqv, env.Q, xlabel="ka", ylabel="TSCS", legend=false),
            layout=(1, 2),
            legend=false
            )

        frame(a, p)
        env(policy(env))
    end

    gif(a, path, fps=fps)
end
