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
