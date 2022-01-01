using DesignEnvironments
using ReinforcementLearning
using Plots
using ProgressMeter

M = 14

design = CoreConfiguration(M = M, plane_size = 15.0, vel_decay=1.0)

a = Animation()

action = randn(2* M) * 3.0
design(action)
action = zeros(2* M)

for i in 1:100
    frame(a, img(design))
    design(action)
end

gif(a, "anim.mp4", fps=20)
closeall()
# savefig(img(design), "config.png")