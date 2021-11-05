using DesignEnvironments
using ReinforcementLearning
# ENV["GKSwstype"] = "nul"
import Flux

env = CylinderEnv(M=10, continuous=true, grid_size=10.0)

traj = CircularArraySARTTrajectory(
    capacity = 1000,
    # state = Matrix{Float64} => (env.params.M, 2),
    state = Tuple{Matrix{Float64}, Vector{Float64}} => (2,),
    action = Vector{Float64} => size(action_space(env))
    )

sampler = NStepBatchSampler{SARTS}(;
        γ = 0.9,
        n = 1,
        stack_size = nothing,
        batch_size = 10,
    )
# cat(coords..., dims=3)

policy = RandomPolicy(action_space(env))

reset!(env)
while !is_terminated(env)

    action = policy(env)
    state = (get_coords(env), env.Q)
    push!(traj; state=state, action=action)
    env(action)

    state = (get_coords(env), env.Q)
    push!(
        traj;
        reward=reward(env),
        state=state,
        terminal=is_terminated(env)
        )
end
