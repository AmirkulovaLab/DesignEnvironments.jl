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

    a = Animation()

    while !is_terminated(env)
        env(policy(env))
        frame(a, img(env))
    end

    gif(a, path, fps=fps)
end
