module DesignEnvironments

export DE
const DE = DesignEnvironments

using ReinforcementLearning
using IntervalSets
using SpecialFunctions
using BlockDiagonals: BlockDiagonal
using LinearAlgebra
using Plots
using ProgressMeter
using Distributions
include("envs/envs.jl")

end
