module DesignEnvironments

export DE
const DE = DesignEnvironments

using ReinforcementLearning
using IntervalSets: ClosedInterval
using SpecialFunctions: besselh, besselj
using BlockDiagonals: BlockDiagonal
using LinearAlgebra
using Plots
using Distributions

include("core/all.jl")
include("designs/all.jl")
include("objectives/all.jl")
include("common/all.jl")

end
