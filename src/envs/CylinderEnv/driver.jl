include("experemental.jl")

M = 10

env = CylinderEnv(M=M)

action = env |> action_space |> rand
env(action)
