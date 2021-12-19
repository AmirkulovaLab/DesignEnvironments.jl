export CylinderEnv, action_space, state_space, reset!, state

mutable struct CylinderEnv <: AbstractEnv
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
        max_vel = 0.2,
        vel_decay = 0.7,
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
        max_vel = max_vel,
        vel_decay = vel_decay,
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

function RLBase.reward(env::CylinderEnv)
    return - env.Q_RMS + env.penalty * env.collision_penalty
end

function RLBase.action_space(env::CylinderEnv)
    if env.continuous
        ## in the case of continuous actions the actino space will be a vector
        step_size = env.config.max_vel
        action_space = Space([-step_size..step_size for _ in 1:(2 * env.config.M)])
    else
        ## in the case of discrete actions it will be an integer within a range.
        action_space = Base.OneTo(4 * env.config.M)
    end
end

function RLBase.state_space(env::CylinderEnv)
    return Space([-Inf..Inf for _ in Base.OneTo(6 * env.config.M + env.nfreq)])
end

function RLBase.reset!(env::CylinderEnv)
    reset_config!(env.config)
    objective(env)
    env.timestep = 0
end

function RLBase.is_terminated(env::CylinderEnv)
    return env.timestep >= env.episode_length
end

coords_to_x(coords::Matrix{Float64})::Vector{Float64} = reshape(coords', length(coords))
x_to_coords(x::AbstractArray{Float64})::Matrix{Float64} = reshape(x, 2, Int(length(x)/2))'

function RLBase.state(env::CylinderEnv)
    pos = coords_to_x(env.config.pos)
    vel = coords_to_x(env.config.vel)
    return vcat(pos, vel, env.qV, env.Q)
end

"""
Converts an integer action into an action matrix which can be
added to the configuration.
"""
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

function objective(env::CylinderEnv)
    x = coords_to_x(env.config.pos)
    env.Q_RMS, env.qV, env.Q = TSCS(x, env.k0amax, env.k0amin, env.nfreq)
end

function (env::CylinderEnv)(action)
    env.timestep += 1

    if !env.continuous
        ## convert discrete action into vector
        action = continuous_action(env, action)
    else
        action = x_to_coords(action)
    end

    env.penalty = 0.0

    env.config(action)
    collisions = get_collisions(env.config)

    env.penalty = - (size(collisions[1], 1) + sum(collisions[2]))
    objective(env)
end

function render(env::CylinderEnv, policy::AbstractPolicy, name::String)
    reset!(env)

    freqv = range(env.k0amin, env.k0amax, length=env.nfreq) |> collect

    a = Animation()
    prog = ProgressUnknown("Working hard:", spinner=true)

    initial = Dict(:rms=>env.Q_RMS, :coords=>env.config.pos)
    optimal = deepcopy(initial)
    max_tscs = maximum(env.Q)

    while !is_terminated(env)
        ProgressMeter.next!(prog)

        ## this plot is the image which displays the current state of the environment
        plot_1 = img(env.config)

        ## this plot shows the scattering pattern (TSCS) produced by the current configuration
        plot_2 = plot(
            freqv, env.Q,
            xlabel="ka", ylabel="TSCS",
            xlim=(freqv[1], freqv[end]),
            ylim=(0, max_tscs), legend=false)

        ## create a side by side plot containing two subplots
        big_plot = plot(
            plot_1,
            plot_2,
            layout=@layout([a{0.6w} b]),
            )

        ## ust this plot as a frame in the animation
        frame(a, big_plot)

        ## apply the policy's action to the environment
        action = policy(env)
        env(action)

        if env.Q_RMS < optimal[:rms]
            optimal[:rms] = env.Q_RMS
            optimal[:coords] = env.config.pos
        end
    end
    ProgressMeter.finish!(prog)

    ## convert collections of images into gif
    gif(a, name * ".mp4", fps=20)
    closeall()

    return (initial, optimal)
end