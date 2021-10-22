module DesignEnvironments

export CylinderEnv

using ReinforcementLearning
using IntervalSets
using SpecialFunctions
using BlockDiagonals
using LinearAlgebra

include("objective.jl")
include("env.jl")

end
