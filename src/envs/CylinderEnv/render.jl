export img, render

#=
    generates a cylindrical shape at specific coordinates with given radius
=#
function cylinder(x, y, r=1; n=30)
    θ = 0:360÷n:360
    Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
end

function img(env::CylinderEnv)

    coords = get_coords(env)
    x, y = coords[:, 1], coords[:, 2]

    ## create a vector of cylinder objects
    cylinders = cylinder.(x, y)

    ## get bounds of grid
    bounds = env.params.grid_size + 1

    p = plot(
        cylinders;
        aspect_ratio=:equal,
        legend=false,
        color=:black,
        xlim=(-bounds, bounds),
        ylim=(-bounds, bounds)
        )

    return p
end

function render(
        env::CylinderEnv, policy::AbstractPolicy;
        max_tscs::Float64=5.0, path::String="anim.mp4",
        fps::Int=30
        )

    reset!(env)

    params = get_params(env)

    k0amin = params.k0amin
    k0amax = params.k0amax
    nfreq = params.nfreq

    dif = (k0amax - k0amin) / nfreq
    freqv = range(k0amin, k0amax, length=nfreq) |> collect

    a = Animation()
    prog = ProgressUnknown("Generating animation:", spinner=true)

    while !is_terminated(env)
        ProgressMeter.next!(prog)

        ## this plot is the image which displays the current state of the environment
        plot_1 = img(env)

        ## this plot shows the scattering pattern (TSCS) produced by the current configuration
        plot_2 = plot(
            freqv, get_Q(env),
            xlabel="ka", ylabel="TSCS",
            xlim=(freqv[1], freqv[end]),
            ylim=(0, max_tscs))

        ## create a side by side plot containing two subplots
        big_plot = plot(
            plot_1,
            plot_2,
            layout=@layout([a{0.6w} b]),
            )

        ## ust this plot as a frame in the animation
        frame(a, big_plot)

        ## apply the policy's action to the environment
        env(policy(env))
    end
    ProgressMeter.finish!(prog)

    ## convert collections of images into gif
    gif(a, path, fps=fps)
    closeall()
end
