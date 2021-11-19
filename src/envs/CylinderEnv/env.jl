export CylinderEnv, action_space, state_space, img

mutable struct CylinderEnv <: DesignEnv
    ## static params
    M::Int
    k0amax::Float64
    k0amin::Float64
    nfreq::Int
    episode_length::Int
    step_size::Float64
    grid_size::Float64
    min_distance::Float64
    velocity_decay::Float64
    continuous::Bool
    physics::Bool

    ## spaces
    action_space::Union{Space, Base.OneTo{Int64}}
    state_space::Space

    ## misc
    timestep::Int

    ## placeholder values for design params
    coords::Matrix{Float64}
    velocity::Matrix{Float64}
    radii::Vector{Float64}
    Q_RMS::Float64
    qV::Vector{Float64}
    Q::Vector{Float64}
    penalty::Float64
end

function CylinderEnv(;
        M::Int = 1,
        radius::Float64 = 1.0,
        k0amax::Float64 = 0.5,
        k0amin::Float64 = 0.3,
        nfreq::Int = 10,
        episode_length::Int = 100,
        step_size::Float64 = 0.5,
        grid_size::Float64 = 5.0,
        min_distance::Float64 = 0.1,
        velocity_decay::Float64=0.9,
        continuous::Bool = true,
        physics::Bool = true
        )

    ## the dimention of the design vector (x) is double the number of cylinders
    x_dim = 2 * M

    if continuous
        ## in the case of continuous actions the actino space will be a vector
        action_space = Space([-step_size..step_size for _ in Base.OneTo(x_dim)])
    else
        ## in the case of discrete actions it will be an integer within a range.
        action_space = Base.OneTo(4 * M)
    end

    #=
        assuming the state of the environment will be represented as a vector the
        content of that vector will be:
            - design vector (x)
            - gradient of design vector (x)
            - Q vector (TSCS at each ka)
            - single float representing the current timestep
    =#

    ## add extra dimentions for velocity vector if physics
    state_dim = (2 + physics) * x_dim + nfreq + 1
    state_space = Space([-Inf..Inf for _ in Base.OneTo(state_dim)])

    timestep = 0
    coords = zeros(M, 2)
    velocity = zeros(M, 2)
    radii = ones(M) .* radius
    Q_RMS = 0.0
    qV = zeros(x_dim)
    Q = zeros(nfreq)
    penalty = 0.0

    ## creating env
    env = CylinderEnv(
        M, k0amax, k0amin, nfreq, episode_length,
        step_size, grid_size, min_distance,
        velocity_decay, continuous, physics,
        action_space, state_space, timestep,
        coords, velocity, radii, Q_RMS, qV, Q, penalty)

    reset!(env)
    return env
end


function reward(env::CylinderEnv)
    return - env.Q_RMS + env.penalty
end

function is_terminated(env::CylinderEnv)
    return env.timestep == env.episode_length
end

function state(env::CylinderEnv)
    coords = coords_to_x(env.coords)

    if env.physics
        velocity = coords_to_x(env.velocity)
        s = vcat(coords, velocity)
    else
        s = coords
    end

    return vcat(s, env.qV, env.Q, env.timestep/env.episode_length)
end

function reset!(env::CylinderEnv)
    env.timestep = 0

    #=
        generates a random uniform configuration of design parameters within the grid
        until a valid one is encountered.
    =#
    while true
        env.coords = (2 * env.grid_size) .* (rand(Float64, env.M, 2) .- 0.5)

        !has_valid_coords(env) || break
    end

    env.velocity = zeros(size(env.coords))

    calculate_objective(env)
end

coords_to_x(coords::Matrix{Float64})::Vector{Float64} = reshape(coords', length(coords))
x_to_coords(x::Vector{Float64})::Matrix{Float64} = reshape(x, 2, Int(length(x)/2))'
#=
    calls the objective function on the current configuration
=#
function calculate_objective(env::CylinderEnv)
    x = coords_to_x(env.coords)
    env.Q_RMS, env.qV, env.Q = TSCS(x, env.k0amax, env.k0amin, env.nfreq)
end

#=
    the purpose of this function is to convert a discrete integer action into
    its equivelant continuous action. each discrete action is capable of moving
    one cylinder in one cardinal direction by a fixed amount.
=#
function continuous_action(env::CylinderEnv, action::Int)
    ## decrement action so that action numbering starts at 0 (instead of 1)
    action -= 1

    action_matrix = zeros(env.M, 2)
    ## getting the cylinder that the current action is adjusting
    cyl = Int(floor(action / 4)) + 1
    ## finding the direction in which that adjustment is being made
    direction = action % 4
    ## finding out if the adjustment is to the x or y axis of cylinder
    axis = (direction % 2) + 1
    ## finding out which direction on the given axis
    sign = Int(floor(direction / 2))
    ## setting the appropriate cylinder and axis equal to the adjustment
    action_matrix[cyl, axis] = (-1)^sign * env.step_size
    ## flattening the action into a vector
    # return coords_to_x(action_matrix)
    return action_matrix
end

#=
    generates a cylindrical shape at specific coordinates with given radius
=#
function cylinder(x, y, r=1; n=30)
    θ = 0:360÷n:360
    Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
end

function get_collisions(env::CylinderEnv)
    ## cartisian product to compare each cylinder to one another
    M = env.M

    idx = 1:M |> collect
    i_idx, j_idx = repeat(idx, outer=M), repeat(idx, inner=M)

    ## outer and inner loop index
    coords = env.coords
    i_coords, j_coords = coords[i_idx, :], coords[j_idx, :]

    ## splitting into coordinates
    x1, y1, x2, y2 = i_coords[:, 1], i_coords[:, 2], j_coords[:, 1], j_coords[:, 2]

    ## calculating distance between each cylinder
    distances = sqrt.((x2 - x1).^2 + (y2 - y1).^2)

    ## if two cylinders are closer than the min distance it is considered a collision
    is_collision = (distances .<= (2 + env.min_distance))

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

function img(env::CylinderEnv)

    coords = env.coords
    x, y = coords[:, 1], coords[:, 2]

    ## create a vector of cylinder objects
    cylinders = cylinder.(x, y)

    p = plot(
        cylinders;
        aspect_ratio=:equal,
        legend=false,
        color=:black,
        xlim=(-env.grid_size, env.grid_size),
        ylim=(-env.grid_size, env.grid_size)
        )

    return p
end

#=
    determines if the current configuration of design parameters is valid.
=#
function has_valid_coords(env::CylinderEnv)::Bool
    within_bounds = false
    overlap = false

    coords = env.coords

    ## first check if all coordinates (x and y) are within the bounds of the grid
    ## radius assumed at 1
    if all(abs.(coords) .< env.grid_size - 1)
        within_bounds = true
        if size(get_collisions(env), 1) > 0
            overlap = true
        end
    end

    return within_bounds & !overlap
end

#=
    defines the effect that the action will have on the environment
=#
function (env::CylinderEnv)(action)
    env.timestep += 1

    ## obtain a deep reference to the current configuration of cylinders
    prev_coords = deepcopy(env.coords)

    if !env.continuous
        ## convert discrete action into vector
        action = continuous_action(env, action)
    else
        action = x_to_coords(action)
    end

    env.penalty = 0.0

    if env.physics
        env.velocity *= env.velocity_decay
        env.velocity += action
        clamp!(env.velocity, -env.step_size, env.step_size)
        env.coords += env.velocity

        wall_collisions = resolve_wall_collisions(env)

        collisions = get_collisions(env)
        resolve_cylinder_collisions(env, collisions)

        env.penalty = - (size(collisions, 1) + sum(wall_collisions))
    else
        ## applying action to configuration
        env.coords += action

        ## check if new configuration is valid
        if !has_valid_coords(env)
            ## setting the current configuration to the last valid one
            env.coords = prev_coords
            ## we want to penalize illegal actions
            env.penalty = -1.0
        end
    end


    ## calculate scattering
    calculate_objective(env)
end

function resolve_wall_collisions(env::CylinderEnv)
    center_bound = env.grid_size .- env.radii

    left_collision = env.coords[:, 1] .< - center_bound
    right_collision = env.coords[:, 1] .> center_bound

    env.coords[left_collision, 1] .= - center_bound[left_collision]
    env.coords[right_collision, 1] .= center_bound[right_collision]
    env.velocity[left_collision .| right_collision, 1] *= - 1

    top_collision = env.coords[:, 2] .> center_bound
    bottom_collision = env.coords[:, 2] .< - center_bound

    env.coords[top_collision, 2] .= center_bound[top_collision]
    env.coords[bottom_collision, 2] .= - center_bound[bottom_collision]
    env.velocity[top_collision .| bottom_collision, 2] *= -1

    return left_collision .+ right_collision .+ top_collision .+ bottom_collision
end


function resolve_cylinder_collisions(env::CylinderEnv, collisions::Matrix)
    ## getting position and velocity of colliding cylinders
    i_pos = env.coords[collisions[:, 1], :]
    j_pos = env.coords[collisions[:, 2], :]
    i_vel = env.velocity[collisions[:, 1], :]
    j_vel = env.velocity[collisions[:, 2], :]

    ## updating velocity of cylinder i
    C_pos = i_pos - j_pos
    C_vel = i_vel - j_vel
    inner_prod = sum(C_vel .* C_pos, dims=2)
    C_pos_norm = C_pos[:, 1].^2 + C_pos[:, 2].^2
    coef = inner_prod ./ C_pos_norm
    env.velocity[collisions[:, 1], :] .-= coef .* C_pos

    ## updating velocity of cylinder j
    C_pos = j_pos - i_pos
    C_vel = j_vel - i_vel
    inner_prod = sum(C_vel .* C_pos, dims=2)
    C_pos_norm = C_pos[:, 1].^2 + C_pos[:, 2].^2
    coef = inner_prod ./ C_pos_norm
    env.velocity[collisions[:, 2], :] .-= coef .* C_pos

    ## Preventing overlaps
    angle = atan.(C_pos[:, 2], C_pos[:, 1])
    midpoint = i_pos + 0.5 * C_pos

    boundry = (2 .* env.radii .+ env.min_distance) ./ 2

    ## Updating x and y coordinates of cylinder i
    env.coords[collisions[:, 1], 1] = midpoint[:, 1] - (cos.(angle) .* boundry[collisions[:, 1]])
    env.coords[collisions[:, 1], 2] = midpoint[:, 2] - (sin.(angle) .* boundry[collisions[:, 1]])

    ## Updating x and y coordinates of cylinder j
    env.coords[collisions[:, 2], 1] = midpoint[:, 1] + (cos.(angle) .* boundry[collisions[:, 2]])
    env.coords[collisions[:, 2], 2] = midpoint[:, 2] + (sin.(angle) .* boundry[collisions[:, 2]])
end
