export CoreConfiguration

function bulid_core(M::Int)::Matrix
    core = [[0.0, 0.0]]
    radius = 2.1
    theta = 2*pi/M

    for i in 1:M
        pos = [radius * cos(theta * i), radius * sin(theta * i)]
        push!(core, pos)
    end

    return hcat(core...)'
end

function hex_core(
        plane::Disc, 
        max_vel::Float64, 
        vel_decay::Float64, 
        min_distance::Float64)

    core = bulid_core(6)

    config = Configuration(
        core,
        zeros(size(core)),
        ones(size(core, 1)),
        plane,
        max_vel,
        vel_decay,
        min_distance)
    
    return config
end

mutable struct CoreConfiguration <: AbstractDesign
    core::Configuration
    config::Configuration
end

function CoreConfiguration(;
    M::Int, plane::Disc,
    max_vel::Float64=MAX_VEL, 
    vel_decay::Float64=VEL_DECAY, 
    min_distance::Float64=MIN_DISTANCE)

    core = hex_core(plane, max_vel, vel_decay, min_distance)
    plane.inner_radius = maximum(core.pos[:, 1]) + core.radii[1] + min_distance

    config = Configuration(
        M = M,
        plane = plane,
        max_vel = max_vel,
        vel_decay = vel_decay,
        min_distance = min_distance)
    
    design = CoreConfiguration(core, config)
    
    reset_design!(design)
    return design
end

"""
Reseting the coordinates of the Configuration of cylinders to be in random positions
around the core.
"""
function reset_design!(design::CoreConfiguration)
    reset_design!(design.config)
end