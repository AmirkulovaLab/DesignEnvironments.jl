"""
A cylinder which is located in a two dimentional plane.
The cylinder has a position, velocity and radius.

Cylinder is randomly initialized anywhere in the plane.

## Example
```
grid_size = 6.0
radius = 1.0

Cylinder(grid_size, radius)
```
"""
mutable struct Cylinder
    position::Vector{Float64}
    velocity::Vector{Float64}
    radius::Float64
end

function Cylinder(grid_size::Float64, radius::Float64)
    position = (2 * grid_size) .* (rand(Float64, 2) .- 0.5)
    velocity = zero(position)
    return Cylinder(position, velocity, radius)
end

function polar_coordinates(cyl::Cylinder; n::Int=30)
    θ = 0:360÷n:360
    return cyl.radius*sind.(θ) .+ cyl.position[1], cyl.radius*cosd.(θ) .+ cyl.position[2]
end

"""
Generate a configuration of randomly positioned M cylinders with given radii / radius

## Example
```
config(10, 6.0, 1.0)
```
"""
function config(M::Int, grid_size::Float64, radii::Vector{Float64})
    config = Vector{Cylinder}(undef, M)

    for i in 1:length(config)
        config[i] = Cylinder(grid_size, radii[i])
    end

    return config
end

function config(M::Int, grid_size::Float64, radius::Float64)
    return config(M, grid_size, fill(radius, M))
end

"""
Gathers the coordinates of each cylinder into a matrix of size:
(M, 2)
"""
function get_positions(cylinders::Vector{Cylinder})::Matrix{Float64}
    positions = getproperty.(cylinders, :position)
    return reduce(hcat, positions)'
end

function set_positions!(cylinders::Vector{Cylinder}, coords::Matrix{Float64})
    coords = [coords[i, :] for i in 1:size(coords, 1)]
    setproperty!.(cylinders, :position, coords)
end

"""
Gathers the velocities of each cylinder into a matrix of size:
(M, 2)
"""
function get_velocities(cylinders::Vector{Cylinder})::Matrix{Float64}
    velocities = getproperty.(cylinders, :velocity)
    return reduce(hcat, velocities)'
end

function set_velocities!(cylinders::Vector{Cylinder}, vel::Matrix{Float64})
    vel = [vel[i, :] for i in 1:size(vel, 1)]
    setproperty!.(cylinders, :velocity, vel)
end

"""
Gathers the radii of all cylinders into a vector of size:
(M,)
"""
function get_radii(cylinders::Vector{Cylinder})
    return getproperty.(cylinders, :radius)
end

function matrix_to_vector(m::Matrix{Float64})::Vector{Float64}
    return reshape(m', length(m))
end

function vector_to_matrix(v::AbstractArray{Float64})::Matrix{Float64}
    return reshape(v, 2, Int(length(v)/2))'
end

"""
Computes the collisions produced by a configuration of cylinders with radii

## Example
```
collisions(coords, radii, min_distance)
```

```
6×2 Matrix{Int64}:
 1   9
 1  10
 2  10
 3   5
 4   7
 9  10
```
"""
function get_collisions(cylinders::Vector{Cylinder}, min_distance::Float64)
    coords = get_positions(cylinders)
    radii = get_radii(cylinders)

    ## cartisian product to compare each cylinder to one another
    M = size(coords, 1)

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
    threshhold_distances = i_radii .+ j_radii .+ min_distance

    ## if two cylinders are closer than the min distance it is considered a collision
    is_collision = (distances .<= threshhold_distances)

    ## we set collision equal to false for comparisons of a cylinder to itself
    same_cylinder_idx = 1:M+1:M^2

    is_collision[same_cylinder_idx] .= false

    ## cylinder index for the ith and jth colliding cylinders
    i_cyl, j_cyl = i_idx[is_collision], j_idx[is_collision]

    ## zipping colliding cylinders together [(a, b), (c, d), (e, f), ...]
    collisions = hcat(i_cyl, j_cyl)
    collisions = sort(collisions, dims=2)
    collisions = unique(collisions, dims=1)

    return collisions
end

"""
Determines if the current configuration of design parameters is valid.
"""
function valid_coords(cylinders::Vector{Cylinder}, grid_size::Float64, min_distance::Float64)
    within_bounds = false
    overlap = false
    coords = get_positions(cylinders)
    radii = get_radii(cylinders)

    ## first check if all coordinates (x and y) are within the bounds of the grid
    if all(abs.(coords) .< grid_size .- radii)
        within_bounds = true
        collisions = get_collisions(cylinders, min_distance)

        if size(collisions, 1) > 0
            overlap = true
        end
    end

    return within_bounds & !overlap
end

"""
Generates a valid configuration of M cylinders with given radii, satisfying
the requirements of grid_size and min_distance.

## Example
```
valid_config(M, grid_size, radius, min_distance)
```
"""
function valid_config(M::Int, grid_size::Float64, r::Union{Vector{Float64}, Float64}, min_distance::Float64)
    cylinders = Vector{Cylinder}(undef, M)

    while true
        cylinders = config(M, grid_size, r)
        !valid_coords(cylinders, grid_size, min_distance) || break
    end

    return cylinders
end


function resolve_cylinder_collisions!(cylinders::Vector{Cylinder}, collisions::Matrix{Int}, min_distance::Float64)
    coords = get_positions(cylinders)
    vel = get_velocities(cylinders)
    radii = get_radii(cylinders)

    ## getting position and velocity of colliding cylinders
    i_pos = coords[collisions[:, 1], :]
    j_pos = coords[collisions[:, 2], :]
    i_vel = vel[collisions[:, 1], :]
    j_vel = vel[collisions[:, 2], :]

    ## updating velocity of cylinder i
    C_pos = i_pos - j_pos
    C_vel = i_vel - j_vel
    inner_prod = sum(C_vel .* C_pos, dims=2)
    C_pos_norm = C_pos[:, 1].^2 + C_pos[:, 2].^2
    coef = inner_prod ./ C_pos_norm
    new_i_vel = i_vel .- coef .* C_pos
    set_velocities!(cylinders[collisions[:, 1]], new_i_vel)
    # new_i_vel = [new_i_vel[i, :] for i in 1:size(new_i_vel, 1)]
    # setproperty!.(cylinders[collisions[:, 1]], :velocity, new_i_vel)

    ## updating velocity of cylinder j
    C_pos = j_pos - i_pos
    C_vel = j_vel - i_vel
    inner_prod = sum(C_vel .* C_pos, dims=2)
    C_pos_norm = C_pos[:, 1].^2 + C_pos[:, 2].^2
    coef = inner_prod ./ C_pos_norm
    new_j_vel = j_vel .- coef .* C_pos
    set_velocities!(cylinders[collisions[:, 2]], new_j_vel)
    # new_j_vel = [new_j_vel[i, :] for i in 1:size(new_j_vel, 1)]
    # setproperty!.(cylinders[collisions[:, 2]], :velocity, new_j_vel)

    ## Preventing overlaps
    angle = atan.(C_pos[:, 2], C_pos[:, 1])
    midpoint = i_pos + 0.5 * C_pos

    boundry = (2 .* radii .+ min_distance) ./ 2

    ## Updating x and y coordinates of cylinder i
    new_xi_coords = midpoint[:, 1] - (cos.(angle) .* boundry[collisions[:, 1]])
    new_yi_coords = midpoint[:, 2] - (sin.(angle) .* boundry[collisions[:, 1]])
    new_i_coords = hcat(new_xi_coords, new_yi_coords)
    set_positions!(cylinders[collisions[:, 1]], new_i_coords)
    # new_i_coords = [new_i_coords[i, :] for i in 1:size(new_i_coords, 1)]
    # setproperty!.(cylinders[collisions[:, 1]], :position, new_i_coords)

    new_xj_coords = midpoint[:, 1] + (cos.(angle) .* boundry[collisions[:, 2]])
    new_yj_coords = midpoint[:, 2] + (sin.(angle) .* boundry[collisions[:, 2]])
    new_j_coords = hcat(new_xj_coords, new_yj_coords)
    set_positions!(cylinders[collisions[:, 2]], new_j_coords)
    # new_j_coords = [new_j_coords[i, :] for i in 1:size(new_j_coords, 1)]
    # setproperty!.(cylinders[collisions[:, 2]], :position, new_j_coords)
end
#
# function resolve_cylinder_collisions(coords::Matrix{Float64}, vel::Matrix{Float64}, radii::Vector{Float64}, collisions::Matrix{Int}, min_distance::Float64)
#     ## getting position and velocity of colliding cylinders
#     i_pos = coords[collisions[:, 1], :]
#     j_pos = coords[collisions[:, 2], :]
#     i_vel = vel[collisions[:, 1], :]
#     j_vel = vel[collisions[:, 2], :]
#
#     ## updating velocity of cylinder i
#     C_pos = i_pos - j_pos
#     C_vel = i_vel - j_vel
#     inner_prod = sum(C_vel .* C_pos, dims=2)
#     C_pos_norm = C_pos[:, 1].^2 + C_pos[:, 2].^2
#     coef = inner_prod ./ C_pos_norm
#     vel[collisions[:, 1], :] .-= coef .* C_pos
#
#     ## updating velocity of cylinder j
#     C_pos = j_pos - i_pos
#     C_vel = j_vel - i_vel
#     inner_prod = sum(C_vel .* C_pos, dims=2)
#     C_pos_norm = C_pos[:, 1].^2 + C_pos[:, 2].^2
#     coef = inner_prod ./ C_pos_norm
#     vel[collisions[:, 2], :] .-= coef .* C_pos
#
#     ## Preventing overlaps
#     angle = atan.(C_pos[:, 2], C_pos[:, 1])
#     midpoint = i_pos + 0.5 * C_pos
#
#     boundry = (2 .* radii .+ min_distance) ./ 2
#
#     ## Updating x and y coordinates of cylinder i
#     coords[collisions[:, 1], 1] = midpoint[:, 1] - (cos.(angle) .* boundry[collisions[:, 1]])
#     coords[collisions[:, 1], 2] = midpoint[:, 2] - (sin.(angle) .* boundry[collisions[:, 1]])
#
#     ## Updating x and y coordinates of cylinder j
#     coords[collisions[:, 2], 1] = midpoint[:, 1] + (cos.(angle) .* boundry[collisions[:, 2]])
#     coords[collisions[:, 2], 2] = midpoint[:, 2] + (sin.(angle) .* boundry[collisions[:, 2]])
# end

function resolve_wall_collisions(cylinders::Vector{Cylinder}, grid_size::Float64)
    coords, vel, radii = get_positions(cylinders), get_velocities(cylinders), get_radii(cylinders)

    ## Setting the boundry for each cylinder
    center_bound = grid_size .- radii

    ## determining which cylinders are colliding with walls
    left_collision = coords[:, 1] .< - center_bound
    right_collision = coords[:, 1] .> center_bound

    coords[left_collision, 1] .= - center_bound[left_collision]
    coords[right_collision, 1] .= center_bound[right_collision]
    vel[left_collision .| right_collision, 1] *= - 1

    ## determining which cylinders are colliding with walls
    top_collision = coords[:, 2] .> center_bound
    bottom_collision = coords[:, 2] .< - center_bound

    coords[top_collision, 2] .= center_bound[top_collision]
    coords[bottom_collision, 2] .= - center_bound[bottom_collision]
    vel[top_collision .| bottom_collision, 2] *= -1

    set_positions!(cylinders, coords)
    set_velocities!(cylinders, vel)
    return left_collision .+ right_collision .+ top_collision .+ bottom_collision
end
