include("cylinder.jl")

mutable struct CylinderEnv
    ## static params
    M::Int
    default_radius::Float64
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

    ## placeholder values for design params
    cylinders::Vector{Cylinder}
    Q_RMS::Float64
    qV::Vector{Float64}
    Q::Vector{Float64}

    ## misc
    penalty::Float64
    timestep::Int
end

function CylinderEnv(;
    M,
    default_radius = 1.0,
    k0amax = 0.5,
    k0amin = 0.3,
    nfreq = 10,
    episode_length = 100,
    step_size = 0.2,
    grid_size = 6.0,
    min_distance = 0.1,
    velocity_decay = 0.7,
    continuous = false,
    physics = true
    )

    cylinders = Vector{Cylinder}(undef, M)
    Q_RMS = 0.0
    qV = zeros(2 * M)
    Q = zeros(nfreq)
    penalty = 0.0
    timestep = 0

    env = CylinderEnv(
        M, default_radius, k0amax, k0amin, nfreq, episode_length,
        step_size, grid_size, min_distance, velocity_decay,
        continuous,physics, cylinders, Q_RMS,
        qV, Q, penalty, timestep)

    reset!(env)
    return env
end

function action_space(env::CylinderEnv)
    if env.continuous
        ## in the case of continuous actions the actino space will be a vector
        action_space = Space([-env.step_size..env.step_size for _ in Base.OneTo(2 * env.M)])
    else
        ## in the case of discrete actions it will be an integer within a range.
        action_space = Base.OneTo(4 * env.M)
    end
end

function reset!(env::CylinderEnv)
    r = env.default_radius ## change this if experementing with radii

    env.cylinders = valid_config(env.M, env.grid_size, r, env.min_distance)
    set_velocities!(env.cylinders, zeros(env.M, 2))
    env.timestep = 0
end

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

function (env::CylinderEnv)(action)
    env.timestep += 1

    ## obtain a deep reference to the current configuration of cylinders
    prev_cylinders = deepcopy(env.cylinders)

    if !env.continuous
        ## convert discrete action into vector
        action = continuous_action(env, action)
    else
        action = vector_to_matrix(action)
    end

    env.penalty = 0.0

    ## current state of cylinders
    coords = get_positions(env.cylinders)
    vel = get_velocities(env.cylinders)
    radii = get_radii(env.cylinders)

    if env.physics
        vel *= env.velocity_decay
        vel += action
        clamp!(vel, -env.step_size, env.step_size)
        coords += vel

        ## assign changes to cylinders
        set_positions!(env.cylinders, coords)
        set_velocities!(env.cylinders, vel)

        ## resolve physical conflicts
        wall_collisions = resolve_wall_collisions(env.cylinders, env.grid_size)
        collisions = get_collisions(get_positions(env.cylinders), radii, env.min_distance)
        resolve_cylinder_collisions!(env.cylinders, collisions, env.min_distance)

        env.penalty = - (size(collisions, 1) + sum(wall_collisions))
    else
        ## applying action to configuration
        coords += action
        set_positions!(cylinders, coords)

        ## check if new configuration is valid
        if !valid_coords(env.cylinders, env.grid_size, env.min_distance)
            ## setting the current configuration to the last valid one
            env.coords = prev_coords
            ## we want to penalize illegal actions
            env.penalty = -1.0
        end
    end
    ## calculate scattering
    # calculate_objective(env)
end
