# define a function that returns a Plots.Shape
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

function Plots.plot(design::Configuration, objective::PressureAmplitude, objective_scale::Tuple)
    p = plot(design, objective, objective_scale)

    xf = pressure_amplitude.xf
    focal_x, focal_y = xf
    focal_point = DE.cylinder(focal_x, focal_y, 1.0)
    grid_size = design.plane.size

    plot!(
        p[1],
        rectangle(grid_size * 2, grid_size * 2, -grid_size, -grid_size),
        color = :black,
        opacity=.1)

    xlim = (
        min(-grid_size, focal_x),
        max(grid_size, focal_x))

    ylim = (
        min(-grid_size, focal_y),
        max(grid_size, focal_y))
    
    plot!(
        p[1], 
        focal_point,
        color = "#e03d6e",
        alpha = 1.0,
        xlim = xlim,
        ylim = ylim)

    plot!(p, size = (1500, 750))
    
    return p
end
