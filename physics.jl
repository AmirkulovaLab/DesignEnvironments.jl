using DesignEnvironments
using ReinforcementLearning
using LinearAlgebra

mutable struct Physics{E} <: AbstractEnv
    env::E
    vel::Matrix{Float64}

    function Physics(env::AbstractEnv)
        vel = zeros(Float64, env.params.M, 2)
        return new{typeof(env)}(env, vel)
    end
end

DesignEnvironments.get_coords(env::Physics{CylinderEnv}) = get_coords(env.env)
DesignEnvironments.img(env::Physics) = img(env.env)

## override ReinforcementLearning functions
RLBase.action_space(env::Physics) = env.env.action_space
RLBase.state_space(env::Physics) = env.env.state_space
RLBase.reward(env::Physics) = env.env.reward
RLBase.is_terminated(env::Physics) = env.env.done
RLBase.state(env::Physics) = DesignEnvironments.get_state(env.env)
RLBase.reset!(env::Physics) = reset!(env.env)

function get_collisions(env::Physics{CylinderEnv})
    ## cartisian product to compare each cylinder to one another
    M = env.env.params.M

    idx = 1:M |> collect
    i_idx = repeat(idx, outer=M)
    j_idx = repeat(idx, inner=M)

    ## outer and inner loop index
    coords = get_coords(env)
    i_coords = coords[i_idx, :]
    j_coords = coords[j_idx, :]

    ## splitting into coordinates
    x1, y1 = i_coords[:, 1], i_coords[:, 2]
    x2, y2 = j_coords[:, 1], j_coords[:, 2]

    ## Calculating distance between each cylinder
    distances = sqrt.((x2 - x1).^2 + (y2 - y1).^2)

    ## If two cylinders are closer than the min distance it is considered a collision
    is_collision = (distances .<= 2.1)

    ## We set collision equal to false for comparisons of a cylinder to itself
    same_cylinder_idx = 1:M+1:M^2

    is_collision[same_cylinder_idx] .= false

    ## Cylinder index for the ith and jth colliding cylinders
    i_cyl, j_cyl = i_idx[is_collision], j_idx[is_collision]

    ## Zipping colliding cylinders together [(a, b), (c, d), (e, f), ...]
    collisions = hcat(i_cyl, j_cyl)
    collisions = sort(collisions, dims=2)
    collisions = unique(collisions, dims=1)

    return collisions
end

function resolve_cylinder_collisions(env::Physics{CylinderEnv})
    collisions = get_collisions(env)
    coords = get_coords(env)
    ## Getting position and velocity of colliding cylinders
    i_pos = coords[collisions[:, 1], :]
    j_pos = coords[collisions[:, 2], :]
    i_vel = env.vel[collisions[:, 1], :]
    j_vel = env.vel[collisions[:, 2], :]

    ## Updating velocity of cylinder i
    C_pos = i_pos - j_pos
    C_vel = i_vel - j_vel
    inner_prod = sum(C_vel .* C_pos, dims=2)
    C_pos_norm = C_pos[:, 1].^2 + C_pos[:, 2].^2
    coef = inner_prod ./ C_pos_norm

    display(coef)

    # inner_prod = np.einsum('ij,ij->i', i_vel-j_vel, C)
    # coef = np.expand_dims(inner_prod / (np.linalg.norm(C, axis=-1) ** 2), -1)
    # self.vel[collisions[:, 0]] -= coef * C
    # ## Updating velocity of cylinder j
    # C = j_pos - i_pos
    # inner_prod = np.einsum('ij,ij->i', j_vel-i_vel, C)
    # coef = np.expand_dims(inner_prod / (np.linalg.norm(C, axis=-1) ** 2), -1)
    # self.vel[collisions[:, 1]] -= coef * C

end

env = Physics(CylinderEnv(M=10, grid_size=10.0))
env.env.x = (2 * env.env.params.grid_size) .* (rand(Float64, 2 * env.env.params.M) .- 0.5)
resolve_cylinder_collisions(env)
display(img(env))

# env = PhysicsCylinderEnv(
#         CylinderEnv(M=10, grid_size=10.0)
#         )
#
#
# resolve_cylinder_collisions(env)
# display(img(env))
