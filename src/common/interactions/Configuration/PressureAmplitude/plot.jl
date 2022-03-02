# define a function that returns a Plots.Shape
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

function Plots.plot(design::Configuration, objective::PressureAmplitude, objective_scale::Tuple)
    p = plot(
        plot(design), 
        plot(objective, objective_scale), 
        layout=@layout([a{0.6w} b]))

    xf = objective.xf
    focal_x, focal_y = xf
    focal_r = 1.0
    focal_point = DE.cylinder(focal_x, focal_y, focal_r)
    grid_size = design.plane.size

    plot!(
        p[1],
        rectangle(grid_size * 2, grid_size * 2, -grid_size, -grid_size),
        color = :black,
        opacity = .1)

    focal_bound = focal_r * 2
    xlim = (
        min(-grid_size, focal_x - focal_bound),
        max(grid_size, focal_x + focal_bound))

    ylim = (
        min(-grid_size, focal_y - focal_bound),
        max(grid_size, focal_y + focal_bound))

    concentration = min(
        1.0,
        metric(objective) / objective_scale[2])

    plot!(
        p[1], 
        focal_point,
        color = "#e03d6e",
        alpha = concentration,
        xlim = xlim,
        ylim = ylim)

    plot!(p, size = (1500, 750))
    
    return p
end
