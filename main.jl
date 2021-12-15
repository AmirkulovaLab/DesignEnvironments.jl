ENV["GKSwstype"] = "nul"
using DesignEnvironments
using Plots
using ReinforcementLearning
RLBase.reward(env::CylinderEnv) = - env.Q_RMS - env.penalty * env.collision_penalty

function objective(env::CylinderEnv)
    x = coords_to_x(env.config.pos)
    env.Q_RMS, env.qV, env.Q = TSCS(x, env.k0amax, env.k0amin, env.nfreq)
end

env = CylinderEnv(
    M=7, 
    plane=Square(10.0),
    nfreq=30,
    k0amax=3.0)

@time objective(env)

p = plot(env.Q)
savefig(p, "TSCS")
savefig(img(env.config), "config")