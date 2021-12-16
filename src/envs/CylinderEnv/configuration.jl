export Configuration, AbstractPlane, Square, Hexagon, img

"""
Defines a type for the plane on which cylinders exist. This type is used to
dispatch the get_wall_collisions function differently for different shapes:
    - Square
    - Hexagon
    - Circle
"""
abstract type AbstractPlane end

struct Square <: AbstractPlane
    size::Float64
end

struct Hexagon <: AbstractPlane
    size::Float64
end

const MAX_VEL = 0.2
const VEL_DECAY = 0.7
const MIN_DISTANCE = 0.1

"""
Simulates a configuration of cylindrical scatterers on a
two dimentional plane.
"""
mutable struct Configuration{P}
    M::Int
    plane::AbstractPlane
    max_vel::Float64
    vel_decay::Float64
    min_distance::Float64
    pos::Matrix
    vel::Matrix
    radii::Vector
end

function Configuration(
        M::Int, plane::AbstractPlane,
        max_vel::Float64, vel_decay::Float64, min_distance::Float64,
        pos::Matrix, vel::Matrix, radii::Vector)

    config = Configuration{typeof(plane)}(M, plane, max_vel, vel_decay, min_distance, pos, vel, radii)
    return config
end

function Configuration(;
        M::Int, plane::AbstractPlane,
        max_vel::Float64=MAX_VEL, vel_decay::Float64=VEL_DECAY, min_distance::Float64=MIN_DISTANCE)

    ## generate starting positions, velocities, and radii
    pos = zeros(M, 2)
    vel = zeros(M, 2)
    radii = ones(M)
    # radii = rand(1:0.01:2, M)

    config = Configuration(M, plane, max_vel, vel_decay, min_distance, pos, vel, radii)

    return config
end

function Configuration(
        pos::Matrix, vel::Matrix, radii::Vector, plane::AbstractPlane,
        max_vel::Float64, vel_decay::Float64, min_distance::Float64)

    M = size(pos, 1)
    config = Configuration(M, plane, max_vel, vel_decay, min_distance, pos, vel, radii)

    return config
end

"""
Generates a set of M random coordinates on a square plane.
"""
function random_pos(config::Configuration)
    center_bound = config.plane.size .- config.radii
    pos = [rand(Uniform(-center_bound[i], center_bound[i]), 2) for i in 1:length(center_bound)]
    return Matrix(hcat(pos...)')
end

function Base.getindex(config::Configuration, I)
    pos = config.pos[I, :]
    vel = config.vel[I, :]
    radii = config.radii[I]


    return Configuration(
        pos, vel, radii, config.plane,
        config.max_vel, config.vel_decay, config.min_distance)
end

"""
Resets the configuration to an initial random state. Sets vel equal to zero
for all cylinders.
"""
function reset_config!(config::Configuration)
    while true
        config.pos = random_pos(config)
        valid_config = size(get_cylinder_collisions(config), 1) == 0
        !valid_config || break
    end

    config.vel = zeros(config.M, 2)
end

function merge_configs(configs...)
    pos = vcat(getfield.(configs, :pos)...)
    vel = vcat(getfield.(configs, :vel)...)
    radii = vcat(getfield.(configs, :radii)...)
    plane = configs[1].plane
    max_vel = configs[1].max_vel
    vel_decay = configs[1].vel_decay
    min_distance = configs[1].min_distance

    return Configuration(
        pos, vel, radii, plane,
        max_vel, vel_decay, min_distance)
end

"""
Calculates the collisions of overlaping cylinders in the configuration.
"""
function get_cylinder_collisions(config::Configuration)
    ## extract some variables
    M = config.M
    coords = config.pos
    radii = config.radii
    ## get indexes of opposing cylinders involved in collision
    idx = collect(1:M)
    i_idx, j_idx = repeat(idx, outer=M), repeat(idx, inner=M)

    ## outer and inner loop index
    i_coords, j_coords = coords[i_idx, :], coords[j_idx, :]
    i_radii, j_radii = radii[i_idx], radii[j_idx]

    ## splitting into coordinates
    x1, y1, x2, y2 = i_coords[:, 1], i_coords[:, 2], j_coords[:, 1], j_coords[:, 2]

    ## calculating distance between each cylinder
    distances = sqrt.((x2 - x1).^2 + (y2 - y1).^2)

    ## calculating threshhold distances for each cylinder comparison
    threshhold_distances = i_radii .+ j_radii .+ config.min_distance

    ## if two cylinders are closer than the min distance it is considered a collision
    is_collision = (distances .<= threshhold_distances)

    ## we set collision equal to false for comparisons of a cylinder to itself
    same_cylinder_idx = 1:M+1:M^2
    is_collision[same_cylinder_idx] .= false

    ## cylinder index for the ith and jth colliding cylinders
    i_cyl, j_cyl = i_idx[is_collision], j_idx[is_collision]

    ## zipping colliding cylinders together [(a, b), (c, d), (e, f), ...]
    collisions = unique(sort(hcat(i_cyl, j_cyl), dims=2), dims=1)
    return collisions
end

function resolve_cylinder_collisions!(config::Configuration, collisions::Matrix)
    coords = config.pos
    vel = config.vel
    radii = config.radii

    ## getting position and vel of colliding cylinders
    i_pos = coords[collisions[:, 1], :]
    j_pos = coords[collisions[:, 2], :]
    i_vel = vel[collisions[:, 1], :]
    j_vel = vel[collisions[:, 2], :]

    ## updating vel of cylinder i
    C_pos = i_pos - j_pos
    C_vel = i_vel - j_vel
    inner_prod = sum(C_vel .* C_pos, dims=2)
    C_pos_norm = C_pos[:, 1].^2 + C_pos[:, 2].^2
    coef = inner_prod ./ C_pos_norm
    config.vel[collisions[:, 1], :] .-= coef .* C_pos

    ## updating vel of cylinder j
    C_pos = j_pos - i_pos
    C_vel = j_vel - i_vel
    inner_prod = sum(C_vel .* C_pos, dims=2)
    C_pos_norm = C_pos[:, 1].^2 + C_pos[:, 2].^2
    coef = inner_prod ./ C_pos_norm
    config.vel[collisions[:, 2], :] .-= coef .* C_pos

    ## Preventing overlaps
    angle = atan.(C_pos[:, 2], C_pos[:, 1])
    midpoint = i_pos + 0.5 * C_pos

    boundry = (2 .* radii .+ config.min_distance) ./ 2

    ## Updating x and y coordinates of cylinder i
    config.pos[collisions[:, 1], 1] = midpoint[:, 1] - (cos.(angle) .* boundry[collisions[:, 1]])
    config.pos[collisions[:, 1], 2] = midpoint[:, 2] - (sin.(angle) .* boundry[collisions[:, 1]])

    ## Updating x and y coordinates of cylinder j
    config.pos[collisions[:, 2], 1] = midpoint[:, 1] + (cos.(angle) .* boundry[collisions[:, 2]])
    config.pos[collisions[:, 2], 2] = midpoint[:, 2] + (sin.(angle) .* boundry[collisions[:, 2]])
end

function get_wall_collisions(config::Configuration{Square})
    coords = config.pos
    ## Setting the boundry for each cylinder
    center_bound = config.plane.size .- config.radii
    ## determining which cylinders are colliding with walls
    left_collision = coords[:, 1] .< - center_bound
    right_collision = coords[:, 1] .> center_bound
    ## determining which cylinders are colliding with walls
    top_collision = coords[:, 2] .> center_bound
    bottom_collision = coords[:, 2] .< - center_bound

    return hcat(left_collision, right_collision, top_collision, bottom_collision)
end

function resolve_wall_collisions!(config::Configuration{Square}, collisions::BitMatrix)
    radii = config.radii
    grid_size = config.plane.size

    ## Setting the boundry for each cylinder
    center_bound = grid_size .- radii

    ## extracting vectors from collision matrix
    left_collision = collisions[:, 1]
    right_collision = collisions[:, 2]
    top_collision = collisions[:, 3]
    bottom_collision = collisions[:, 4]

    ## updating x position and vel
    config.pos[left_collision, 1] .= - center_bound[left_collision]
    config.pos[right_collision, 1] .= center_bound[right_collision]
    config.vel[left_collision .| right_collision, 1] *= - 1

    ## updating y position and vel
    config.pos[top_collision, 2] .= center_bound[top_collision]
    config.pos[bottom_collision, 2] .= - center_bound[bottom_collision]
    config.vel[top_collision .| bottom_collision, 2] *= -1
end

"""
Will use the calculate the cylinder collisions within the configuration as well
as collisions between the cylinders and the walls. Returns both collision sets
together as a tuple of matrix.
"""
function get_collisions(config::Configuration)
    cylinder_collisions = get_cylinder_collisions(config)
    wall_collisions = get_wall_collisions(config)
    return (cylinder_collisions, wall_collisions)
end

"""
Takes in the configuration and its collisions. Resolves overlaps and wall
collisions.
"""
function resolve_collisions!(config::Configuration, collisions::Tuple)
    cylinder_collisions, wall_collisions = collisions
    resolve_cylinder_collisions!(config, cylinder_collisions)
    resolve_wall_collisions!(config, wall_collisions)
end

function (config::Configuration)(action)
    config.vel *= config.vel_decay
    config.vel += action
    clamp!(config.vel, -config.max_vel, config.max_vel)
    config.pos += config.vel

    collisions = get_collisions(config)
    resolve_collisions!(config, collisions)

    return collisions
end

function cylinder(x, y, r; n=30)
    θ = 0:360÷n:360
    Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
end

function img(config::Configuration{Square})
    coords = config.pos
    radii = config.radii
    grid_size = config.plane.size
    x, y = coords[:, 1], coords[:, 2]

    ## create a vector of cylinder objects
    cylinders = cylinder.(x, y, radii)

    p = plot(
        cylinders;
        aspect_ratio=:equal,
        legend=false,
        color=:black,
        xlim=(-grid_size, grid_size),
        ylim=(-grid_size, grid_size))

    return p
end
