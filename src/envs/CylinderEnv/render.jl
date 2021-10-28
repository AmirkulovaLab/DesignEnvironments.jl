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
        max_tscs::Float64=5.0, path::String="anim.mp4",
        fps::Int=30
        )

    reset!(env)

    k0amin = env.params.k0amin
    k0amax = env.params.k0amax
    nfreq = env.params.nfreq

    dif = (k0amax - k0amin) / nfreq
    freqv = range(k0amin, k0amax, length=nfreq) |> collect

    a = Animation()

    while !is_terminated(env)
        plot_1 = img(env)
        plot_2 = plot(
            freqv, env.Q,
            xlabel="ka", ylabel="TSCS",
            xlim=(freqv[1], freqv[end]),
            ylim=(0, max_tscs))

        big_plot = plot(
            plot_1,
            plot_2,
            layout=@layout([a{0.6w} b]),
            )

        frame(a, big_plot)
        env(policy(env))
    end

    gif(a, path, fps=fps)
    closeall()
end
