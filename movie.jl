using Plots
using DesignEnvironments

function circle(x, y, r=1; n=30)
    θ = 0:360÷n:360
    Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
end

kwargs = (
    aspect_ratio=:equal,
    legend=false,
    color=:black,
    xlim=(-5,5),
    ylim=(-5,5))

env = CylinderEnv(M=10)

coords = get_coords(env)

x, y = coords[:, 1], coords[:, 2]

circles = circle.(x, y)

plot(circles; kwargs...)
