export Configuration, Square, Disc

"""
Defines a type for the plane on which cylinders exist. This type is used to
dispatch the `get_wall_collisions` function differently for different shapes.
"""
abstract type AbstractPlane end

"""
A Square plane

# Example
```
plane = Square(15.0)
```
"""
struct Square <: AbstractPlane
    size::Float64
end

"""
Defines the `size` function for `Square` planes. Returns the plane's size.
"""
Base.size(plane::Square) = plane.size

"""
A `Disc` plane with an inner and outer radius. Configurations of cylinders can
exist between these radii.

# Example
```
plane = Disc(2.0, 15.0)
```
"""
mutable struct Disc <: AbstractPlane
    inner_radius::Float64
    outer_radius::Float64
end

"""
Provides an option of only specifying the outer radius. Assumes the inner radius
is zero.

# Example
```
plane = Disc(15.0)
```
"""
function Disc(outer_radius::Float64)
    return Disc(0.0, outer_radius)
end

"""
Returns the outer radius of the `Disc` plane.
"""
Base.size(plane::Disc) = plane.outer_radius

"""
Simulates a configuration of cylindrical scatterers on a
two dimentional plane.

# Parameters
- `M::Int`: number of cylindrical scatterers
- `plane::AbstractPlane`: 2 dimentional surface on which scatterers exist
- `max_vel::Float64`: maximum allowable movement speed of cylinders
- `vel_decay::Float64`: percentage of velocity which remains after each step (0, 1.0)
- `min_distance::Float64`: minimum allowable distance between cylinders
- `pos::Matrix`: (x, y) coordinates of each cylinder
- `vel::Matrix`: (dx, dy) velocity of each cylinder
- `radii::Vector`: radius of each cylinder
"""
mutable struct Configuration{P <: AbstractPlane} <: AbstractDesign
    M::Int
    plane::P
    max_vel::Float64
    vel_decay::Float64
    min_distance::Float64
    pos::Matrix
    vel::Matrix
    radii::Vector
end

"""
Constructor which initializes a randomly generated Configuration.

# Example
```
config = Configuration(
    M = 7,
    plane = Square(15.0),
    max_vel = 0.1,
    vel_decay = 0.7,
    min_distance = 0.1)
```

The arguments `max_vel`, `vel_decay`, and `min_distance` are optional.
"""
function Configuration(
        M::Int, plane::AbstractPlane,
        max_vel::Float64, vel_decay::Float64, min_distance::Float64)

    ## generate starting positions, velocities, and radii
    pos = zeros(M, 2)
    vel = zeros(M, 2)
    radii = ones(M)

    config = Configuration(M, plane, max_vel, vel_decay, min_distance, pos, vel, radii)

    reset_design!(config)
    
    return config
end

"""
Constructor for instantiating a `Configuration` without specifying plane type.
Plane is assumed to be `Square`.
"""
function Configuration(
        M::Int, plane_size::Float64, 
        max_vel::Float64, vel_decay::Float64, min_distance::Float64)

    return Configuration(M, Square(plane_size), max_vel, vel_decay, min_distance)
end

"""
Constructor for defining a Configuration from a preexisting set of coordinates.

# Example
```
pos = zeros(M, 2)
vel = zeros(M, 2)
radii = ones(M)

config = Configuration(
    pos, vel, radii, Square(15.0), 
    0.1, 0.7, 0.1)
```
"""
function Configuration(
        pos::Matrix, vel::Matrix, radii::Vector, plane::AbstractPlane,
        max_vel::Float64, vel_decay::Float64, min_distance::Float64)

    M = size(pos, 1)
    config = Configuration(M, plane, max_vel, vel_decay, min_distance, pos, vel, radii)

    return config
end

"""
Converts a `Matrix` of coordinates `(M, 2)` to a `Vector` of length `M * 2`.
"""
coords_to_x(coords::Matrix)::Vector = reshape(coords', length(coords))

"""
Converts a `Vector` of coordinates to a `Matrix`
"""
x_to_coords(x::AbstractArray)::Matrix = reshape(x, 2, Int(length(x)/2))'

"""
Generates a set of `M` random coordinates on a `Square` plane. Must be  
`Configuration{Square}` type.
"""
function random_pos(config::Configuration{Square})
    center_bound = config.plane.size .- config.radii
    pos = [rand(Uniform(-center_bound[i], center_bound[i]), 2) for i in 1:length(center_bound)]
    return Matrix(hcat(pos...)')
end

"""
Generates a set of M random coordinates on a `Disc` plane. Must be 
`Configuration{Disc}` type.
"""
function random_pos(config::Configuration{Disc})
    inner_radius = config.plane.inner_radius .+ config.radii

    outer_radius = config.plane.outer_radius .- config.radii

    r = [rand(Uniform(i, o)) for (i, o) in zip(inner_radius, outer_radius)]
    θ = rand(Uniform(0, 2 * pi), config.M)

    x = cos.(θ) .* r
    y = sin.(θ) .* r
    
    return hcat(x, y)
end

"""
Resets the `Configuration` to an initial random state. Sets vel equal to zero
for all cylinders.
"""
function reset_design!(config::Configuration)

    ## generate a random valid configuration
    while true
        config.pos = random_pos(config)
        valid_config = size(get_cylinder_collisions(config), 1) == 0
        !valid_config || break
    end

    ## reset velocity to zeros
    config.vel = zeros(config.M, 2)
end

"""
Obtains the position, velocity, and radius of the `i` th cylinder.
"""
function Base.getindex(config::Configuration, I)
    pos = config.pos[I, :]
    vel = config.vel[I, :]
    radii = config.radii[I]


    return Configuration(
        pos, vel, radii, config.plane,
        config.max_vel, config.vel_decay, config.min_distance)
end

"""
Combines many `Configuration` (s) into one.

# Example
```
config = merge_configs(config1, config2, config3)
```
"""
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
Calculates the collisions of overlaping cylinders in the `Configuration`.
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

"""
Resolves any inter-cylinder collisions that are present in the configuration. Returns the
'Configuration' to a valid state.

# Example
```
collisions = get_cylinder_collisions(config)
resolve_cylinder_collisions!(config, collisions)
```
"""
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

"""
Determines the collisions for a `Configuration` on a `Square` plane. Collisions may occur
against any of the four walls.

Must be called on `Configuration{Square}`
"""
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

"""
Determines the wall collisions for a `Configuration`` which exists on a `Disc` plane.
Collisions may occur against the inner radius or outer radius of the disc.

Must be called on `Configuration{Disc}`
"""
function get_wall_collisions(config::Configuration{Disc})
    coords = config.pos

    ## calculating distance of cylinders from origin
    x, y = coords[:, 1], coords[:, 2]
    d = sqrt.(x.^2 + y.^2)

    ## determining if cylinders are less than the inner radius
    inner_collision = (d .- config.radii) .< config.plane.inner_radius
    outer_collision = (d .+ config.radii) .> config.plane.outer_radius

    return hcat(inner_collision, outer_collision)
end

"""
Resolves cylinder-wall collisions by adjusting the position and velocities of the offending
cylinders accordingly.

# Arguments
- `config::Configuration{Square}`: current `Configuration` on a `Square` plane.
- `collisions::BitMatrix`: cylinders colliding with walls.
"""
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
# Arguments
- `config::Configuration{Disc}`
- `collisions::BitMatrix`
"""
function resolve_wall_collisions!(config::Configuration{Disc}, collisions::BitMatrix)
    ## extracting inner and outer collisions from the BitMatrix
    inner_collision = collisions[:, 1]
    outer_collision = collisions[:, 2]

    ## extracting x and y coordinates of all cylinders
    x, y = config.pos[:, 1], config.pos[:, 2]

    ## calculating the new coordinates for cylinders colliding with the inner radius of the Disc
    # d = sqrt.(x[inner_collision] .^ 2 + y[inner_collision] .^ 2)
    # theta = asin.(y[inner_collision] ./ d)

    theta = atan.(y[inner_collision], x[inner_collision])
    r =  config.plane.inner_radius .+ config.radii[inner_collision]

    ## updating the coordinates of inner collisions
    config.pos[inner_collision, 1] .= r .* cos.(theta)
    config.pos[inner_collision, 2] .= r .* sin.(theta)

    ## calculating the new coordinates for cylinders colliding with the outer radius
    # d = sqrt.(x[outer_collision] .^ 2 + y[outer_collision] .^ 2)
    # theta = asin.(y[outer_collision] ./ d)

    theta = atan.(y[outer_collision], x[outer_collision])
    r = config.plane.outer_radius .- config.radii[outer_collision]

    ## updating the coordinates of inner collisions
    config.pos[outer_collision, 1] .= r .* cos.(theta)
    config.pos[outer_collision, 2] .= r .* sin.(theta)

    config.vel[inner_collision .| outer_collision, 1] *= -1
    config.vel[inner_collision .| outer_collision, 2] *= -1
end

"""
Calculates inter-cylinder collisions as well as collisions between the 
cylinders and the walls. Returns both collision sets together as a `Tuple` of `Matrix`.

# Arguments
- `config::Configuration`

# Returns
- `::Tuple`
"""
function get_collisions(config::Configuration)
    cylinder_collisions = get_cylinder_collisions(config)
    wall_collisions = get_wall_collisions(config)
    return (cylinder_collisions, wall_collisions)
end

"""
Takes in the configuration and its collisions. Resolves overlaps and wall
collisions.

# Arguments
- `config::Configuration`
- `collisions::Tuple`
"""
function resolve_collisions!(config::Configuration, collisions::Tuple)
    cylinder_collisions, wall_collisions = collisions
    resolve_cylinder_collisions!(config, cylinder_collisions)
    resolve_wall_collisions!(config, wall_collisions)
end

"""
Generates points in a circular shape at the given (x, y) location with the given radius.

# Arguments
- `x::Float64`
- `y::Float64`
- `r::Float64`

# Returns
- `::Plots.Shape`
"""
function cylinder(x, y, r; n=30)
    θ = 0:360÷n:360
    Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
end

"""
Plots an image of the current state of the `Configuration`.
"""
function Plots.plot(config::Configuration; color=:teal)
    coords = config.pos
    radii = config.radii
    grid_size = size(config.plane)
    x, y = coords[:, 1], coords[:, 2]

    ## create a vector of cylinder objects
    cylinders = cylinder.(x, y, radii)

    p = plot(
        cylinders;
        aspect_ratio=:equal,
        legend=false,
        color=color,
        xlim=(-grid_size, grid_size),
        ylim=(-grid_size, grid_size))

    return p
end

"""
Converts an integer action into an action matrix which can be
added to the configuration.

# Arguments
- `config::Configuration`
- `action::Int`

# Returns
- `::Matrix`

# Example
```
action = continuous_action(config, 5)
```
"""
function continuous_action(config::Configuration, action::Int)
    ## decrement action so that action numbering starts at 0 (instead of 1)
    action -= 1

    action_matrix = zeros(config.M, 2)
    ## getting the cylinder that the current action is adjusting
    cyl = Int(floor(action / 4)) + 1
    ## finding the direction in which that adjustment is being made
    direction = action % 4
    ## finding out if the adjustment is to the x or y axis of cylinder
    axis = (direction % 2) + 1
    ## finding out which direction on the given axis
    sign = Int(floor(direction / 2))
    ## setting the appropriate cylinder and axis equal to the adjustment
    action_matrix[cyl, axis] = (-1)^sign * config.max_vel
    ## flattening the action into a vector
    return action_matrix
end

"""
Function which handles the application of an action to the design. The action can be
either a matrix of velocity adjustments or an integer which represents a discrete action.

The function returns a penalty from the action which is to be minimized in the reward function.

# Arguments
- `action::Union{Vector, Int}`

# Returns
- `::Int`

# Example
```
action = randn(config.M * 2)
penalty = config(action)
```
"""
function (config::Configuration)(action)

    ## convert integer actions to a matrix of velocity adjustments
    if typeof(action) == Int
        action = continuous_action(config, action)
    else
        action = x_to_coords(action)
    end

    ## decay velocity
    config.vel *= config.vel_decay
    ## apply acceleration (action)
    config.vel += action
    ## enforce speed limit
    clamp!(config.vel, -config.max_vel, config.max_vel)
    ## update position
    config.pos += config.vel
    ## determine 
    collisions = get_collisions(config)
    resolve_collisions!(config, collisions)

    ## sum inter cylinder collisions and wall collisions
    penalty = (size(collisions[1], 1) + sum(collisions[2]))
    return penalty
end