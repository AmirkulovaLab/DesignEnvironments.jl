export
    CylinderEnv, action_space, state_space,
    reward, is_terminated, state, reset!,
    get_coords

#=
    immutable struct which holds the static parameters of the cylinder environment
=#
struct CylinderEnvParams
    M::Int
    k0amax::Real
    k0amin::Real
    nfreq::Int
    episode_length::Int
    step_size::Real
    grid_size::Real
    continuous::Bool
end

#=
    this struct holds the dynamic information which describes the state of the
    environment.
=#
mutable struct CylinderEnv <: AbstractEnv
    ## static params
    params::CylinderEnvParams

    ## spaces
    action_space::Union{Space, Base.OneTo{Int64}}
    state_space::Space

    ## misc
    state::Vector
    reward::Real
    done::Bool
    timestep::Int

    ## placeholder values for design params
    x::Vector{Float64}
    Q_RMS::Real
    qV::Vector{Float64}
    Q::Vector{Float64}
end

function CylinderEnv(;
    M::Int = 1,
    k0amax::Real = 0.5,
    k0amin::Real = 0.3,
    nfreq::Int = 10,
    episode_length::Int = 100,
    step_size::Real = 0.5,
    grid_size::Real = 5.0,
    continuous::Bool = true
    )

    ## setting static params
    params = CylinderEnvParams(
        M, k0amax, k0amin, nfreq, episode_length, step_size, grid_size, continuous
        )

    x_dim = 2 * params.M

    if continuous
        action_space = Space([
            ClosedInterval(-params.step_size, params.step_size) for _ in Base.OneTo(x_dim)
            ])
    else
        action_space = Base.OneTo(4 * params.M)
    end

    state_dim = 2 * x_dim + params.nfreq + 1
    state_space = Space([ClosedInterval(-Inf, Inf) for _ in Base.OneTo(state_dim)])

    ## creating env
    env = CylinderEnv(
        params,
        action_space,
        state_space,
        zeros(state_dim), ## state
        0.0, ## reward
        false, ## done
        0, ## timestep
        zeros(x_dim), ## x (design params)
        0.0, ## Q_RMS
        zeros(x_dim), ## qV (gradient)
        zeros(params.nfreq) ## Q (TSCS across nfreq)
        )

    reset!(env)
    return env
end

## various getters
get_state(env::CylinderEnv)::Vector = vcat(env.x, env.qV, env.Q, env.timestep/env.params.episode_length)
get_reward(env::CylinderEnv)::Real = -env.Q_RMS
get_done(env::CylinderEnv)::Bool = env.timestep == env.params.episode_length

## override ReinforcementLearning functions
RLBase.action_space(env::CylinderEnv) = env.action_space
RLBase.state_space(env::CylinderEnv) = env.state_space
RLBase.reward(env::CylinderEnv) = env.reward
RLBase.is_terminated(env::CylinderEnv) = env.done
RLBase.state(env::CylinderEnv) = get_state(env)

#=
    converts the current vector form of the design parameters into a matrix containing
    the coordinates of the cylinders within the grid
=#
get_coords(env::CylinderEnv)::Matrix = transpose(reshape(env.x, 2, env.params.M))

#=
    determines if the current configuration of design parameters is valid.
=#
function has_valid_x(env::CylinderEnv)::Bool
    within_bounds = false
    overlap = false

    coords = get_coords(env)
    if all(abs.(coords) .< env.params.grid_size)
        within_bounds = true
        for i = 1 : env.params.M
            for j = 1 : env.params.M
                if i != j
                    x1, y1 = coords[i, :]
                    x2, y2 = coords[j, :]
                    d = sqrt((x2 - x1)^2 + (y2 - y1)^2)
                    if d <= 2.1
                        overlap = true
                    end
                end
            end
        end
    end

    return within_bounds & !overlap
end

#=
    generates a random uniform configuration of design parameters within the grid
    until a valid one is encountered.
=#
function generate_valid_x(env::CylinderEnv)
    while true
        env.x = (2 * env.params.grid_size) .* (rand(Float64, 2 * env.params.M) .- 0.5)
        !has_valid_x(env) || break
    end
end

#=
    calls the objective function on the current configuration
=#
function calculate_objective(env::CylinderEnv)
    env.Q_RMS, env.qV, env.Q = TSCS(env.x, env.params.k0amax, env.params.k0amin, env.params.nfreq)
end

#=
    resets environment to random starting design
=#
function RLBase.reset!(env::CylinderEnv)
    env.timestep = 0

    ## generate new design parameters
    generate_valid_x(env)

    ## calculate scattering
    calculate_objective(env)

    ## set env info
    env.state = get_state(env)
    env.reward = get_reward(env)
    env.done = false
end

function continuous_action(env::CylinderEnv, action::Int)
    action -= 1

    action_matrix = zeros(env.params.M, 2)
    cyl = Int(floor(action / 4)) + 1
    direction = action % 4
    axis = (direction % 2) + 1
    sign = Int(floor(direction / 2))

    action_matrix[cyl, axis] = (-1)^sign * env.params.step_size
    return reshape(transpose(action_matrix), length(action_matrix))
end

#=
    defines the effect that the action will have on the environment
=#
function (env::CylinderEnv)(action::Union{AbstractArray, Int})
    env.timestep += 1

    prev_x = deepcopy(env.x)

    if !env.params.continuous
        ## convert discrete action into vector
        action = continuous_action(env, action)
    end

    env.x += action

    ## check if new configuration is valid
    if !has_valid_x(env)
        env.x = prev_x
    end

    ## calculate scattering
    calculate_objective(env)

    ## env info
    env.state = get_state(env)
    env.reward = get_reward(env)
    env.done = get_done(env)
end
