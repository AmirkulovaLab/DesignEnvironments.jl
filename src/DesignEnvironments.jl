module DesignEnvironments

export
    CylinderEnv, action_space, state_space
    reward, is_terminated, state, reset!

using ReinforcementLearning
using IntervalSets
using SpecialFunctions
using BlockDiagonals
using LinearAlgebra

include("objective.jl")
include("env.jl")

end
