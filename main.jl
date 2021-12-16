ENV["GKSwstype"] = "nul"
using DesignEnvironments
using Plots
using ReinforcementLearning

env = CylinderEnv(
    M=10, 
    plane=Square(10.0),
    nfreq=20,
    k0amax=1.0,
    max_vel=0.1,
    vel_decay=1.0)

mutable struct InitialActionPolicy <: AbstractPolicy
    action_space::Space
    action::Vector
end

function InitialActionPolicy(action_space::Space)
    action = rand(action_space)
    return InitialActionPolicy(action_space, action)
end

function (p::InitialActionPolicy)(env::CylinderEnv)
    action = deepcopy(p.action)
    p.action = zeros(size(p.action_space))
    return action
end

policy = RandomPolicy(action_space(env))
# policy = InitialActionPolicy(action_space(env))

display("With $(Threads.nthreads()) threads")

DE.render(env, policy, "animation")