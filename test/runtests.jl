using DesignEnvironments
using Test

@testset "DesignEnvironments" begin
    env = CylinderEnv(M=3)
    display(env)
end
