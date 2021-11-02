using DesignEnvironments
using ReinforcementLearning
# ENV["GKSwstype"] = "nul"
using Flux

struct GraphEmb
    emb::Matrix
end

function GraphEmb(emb_dim::Int, graph_dim::Int)
    emb = randn(graph_dim, emb_dim)
    return GraphEmb(emb)
end

Flux.@functor GraphEmb

function (m::GraphEmb)(x)
    emb = m.emb

    x = batched_mul(x, emb)
    x = sum(x, dims=3)
    x = relu.(x)
    return x
end

emb = GraphEmb(10, 2)
x = randn(20, 2, 5)

display(size(emb(x)))

# display(Threads.nthreads())
#
# env = CylinderEnv(M=10, continuous=true, grid_size=10.0)
#
# min_M = 3
# max_M = 10
#
# envs = Vector{CylinderEnv}(undef, 100)
#
# @time begin
#     Threads.@threads for i = 1 : length(envs)
#         envs[i] = CylinderEnv(M=rand(min_M:max_M), continuous=false, grid_size=10.0)
#     end
#
#     Threads.@threads for i = 1 : length(envs)
#         reset!(envs[i])
#     end
# end
