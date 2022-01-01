using DesignEnvironments
using ReinforcementLearning
using Plots
using ProgressMeter

# env = DesignEnvironment(
#     design = Configuration(M = 12, plane = Square(15.0)), 
#     objective = TSCS(k0amax = 1.0, k0amin = 0.3, nfreq = 20),
#     is_continuous = false
#     )

# policy = RandomPolicy(action_space(env))


# a = Animation()
# prog = ProgressUnknown("Working hard:", spinner=true)

# max_tscs = maximum(env.objective.Q)

# while !is_terminated(env)
#     ProgressMeter.next!(prog)

#     env(policy(env))
#     frame(a, img(env))
# end

# ProgressMeter.finish!(prog)

# ## convert collections of images into gif
# gif(a, "better" * ".mp4", fps=20)
# closeall()

function DE.img(design::CoreConfiguration)
    core_img = img(design.core, color=:plum)

    x, y = design.config.pos[:, 1], design.config.pos[:, 2]
    r = design.config.radii

    ## create a vector of cylinder objects
    cylinders = DE.cylinder.(x, y, r)
    plot!(core_img, cylinders; color=:teal)
    return core_img
end

design = CoreConfiguration(
    M = 15,
    plane = Disc(15.0))

DE.get_wall_collisions(design.config)

savefig(img(design), "config.png")

# coords = design.config.pos
# r1 = design.config.radii
# grid_size = design.config.plane.size
# x1, y1 = coords[:, 1], coords[:, 2]

# x2, y2 = design.core.pos[:, 1], design.core.pos[:, 2]
# r2 = design.core.radii

# ## create a vector of cylinder objects
# cylinders1 = DE.cylinder.(x1, y1, r1)
# cylinders2 = DE.cylinder.(x2, y2, r2)

# p = plot(
#     cylinders1;
#     aspect_ratio=:equal,
#     legend=false,
#     color=:teal,
#     xlim=(-grid_size, grid_size),
#     ylim=(-grid_size, grid_size))

# plot!(p, cylinders2; color=:plum)

# savefig(p, "config.png")

# tscs = TSCS(k0amax = 1.0, k0amin = 0.3, nfreq = 20)