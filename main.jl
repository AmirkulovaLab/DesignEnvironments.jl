ENV["GKSwstype"] = "nul"
using DesignEnvironments
using Plots
using ReinforcementLearning

env = CylinderEnv(
    M=10, 
    plane=Square(10.0),
    nfreq=10,
    k0amax=1.0)


for i in 1:10
    reset!(env)
    @time DE.objective(env)
    # savefig(plot(env.Q), "TSCS.png")
    # savefig(img(env.config), "config.png")
end

display("With $(Threads.nthreads()) threads")