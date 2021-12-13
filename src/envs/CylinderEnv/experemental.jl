include("configuration.jl")

mutable struct CylinderEnv
    config::Configuration
    k0amax::Float64
    k0amin::Float64
    nfreq::Int
    episode_length::Int
    continuous::Bool
    collision_penalty::Float64
    penalty::Float64
    timestep::Int

    ## placeholder values for design params
    Q_RMS::Float64
    qV::Vector{Float64}
    Q::Vector{Float64}
end

function CylinderEnv(;
        M,
        plane,
        max_velocity = 0.2,
        velocity_decay = 0.7,
        min_distance = 0.1,
        k0amax = 0.5,
        k0amin = 0.3,
        nfreq = 10,
        episode_length = 100,
        continuous = true,
        collision_penalty = 0.1)

    config = Configuration(
        M = M,
        plane = plane,
        max_velocity = max_velocity,
        velocity_decay = velocity_decay,
        min_distance = min_distance)

    Q_RMS = 0.0
    qV = zeros(2 * M)
    Q = zeros(nfreq)
    penalty = 0.0
    timestep = 0

    env = CylinderEnv(
        config, k0amax, k0amin, nfreq,
        episode_length, continuous,
        collision_penalty, penalty, timestep,
        Q_RMS, qV, Q)

    reset!(env)
    return env
end

function action_space(env::CylinderEnv)
    if env.continuous
        ## in the case of continuous actions the actino space will be a vector
        action_space = Space([-env.step_size..env.step_size for _ in Base.OneTo(2 * env.M)])
    else
        ## in the case of discrete actions it will be an integer within a range.
        action_space = Space([1:(4 * env.M)])
    end
end

function reset!(env::CylinderEnv)
    reset!(env.config)
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
    return action_matrix
end

function (env::CylinderEnv)(action)
    env.timestep += 1

    if !env.continuous
        ## convert discrete action into vector
        action = continuous_action(env, action)
    else
        action = vector_to_matrix(action)
    end

    env.penalty = 0.0

    env.config(action)
    collisions = get_collisions(env.config)

    env.penalty = - (size(collisions[1], 1) + sum(collisions[2]))
end

ENV["GKSwstype"] = "nul"

config = Configuration(
    M = 10,
    plane = Square(8.0),
    velocity_decay=1.0)

action = randn(config.M, 2)
config(action)
action = zeros(config.M, 2)

# a = Animation()

for i in 1:500
    display(config.vel)
    # frame(a, img(config))
    # action = randn(config.M, 2)
    config(action)
end

# gif(a, "test.mp4", fps=20)
# closeall()
