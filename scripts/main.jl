using DesignEnvironments
using ReinforcementLearning

params = Dict(
    :is_continuous => false,
    :design_params => Dict(
        :M => 7,
        :plane_size => 15.0,
        :vel_decay => 0.9
    ),
    :objective_params => Dict(
        :k0amax => 1.0,
        :k0amin => 0.3,
        :nfreq => 15
    )
)

env = DesignEnvironment(
    design = CoreConfiguration(;params[:design_params]...),
    objective = TSCS(;params[:objective_params]...),
    is_continuous=params[:is_continuous]
    )

policy = RandomPolicy(action_space(env))

run(policy, env, StopWhenDone(), TotalRewardPerEpisode())