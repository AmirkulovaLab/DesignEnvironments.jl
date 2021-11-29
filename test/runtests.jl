using DesignEnvironments
using ReinforcementLearning
using Test
ENV["GKSwstype"] = "nul"

@testset "DesignEnvironments" begin
    env = CylinderEnv(M=10, continuous=false)
    policy = RandomPolicy(action_space(env))

    initial, optimal = render(env, policy, path="test.mp4")
    display(initial)
    display(optimal)
end
