using DesignEnvironments
using ReinforcementLearning
using Plots
using ProgressMeter

design = CoreConfiguration(M = 7, plane_size = 15.0)

savefig(img(design), "config.png")