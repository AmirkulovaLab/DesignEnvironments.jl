# ENV["GKSwstype"] = "nul"

using DesignEnvironments
using ReinforcementLearning

env = Physics(CylinderEnv(M=10, grid_size=10, k0amax=2.0, nfreq=20, continuous=false))

# policy = RandomPolicy(action_space(env))
# render(env, policy, path="thing.mp4", max_tscs=30.0)

collisions = DesEnv.get_collisions(env)
DesEnv.resolve_cylinder_collisions(env)

function resolve_wall_collisions(env::Physics{CylinderEnv})
    params = get_params(env)
    coords = get_coords(env)

    ## radius size is assumed to be 1
    center_bound = params.grid_size - 1

    left_collision = coords[:, 1] .< - center_bound
    right_collision = coords[:, 1] .> center_bound

    coords[left_collision, 1] .= - center_bound
    coords[right_collision, 1] .= center_bound
    env.vel[]

    display(left_collision)
    println()
    display(right_collision)
end

resolve_wall_collisions(env)
display(img(env))

#
# def resolve_wall_collisions(self):
#         '''
#         Resolves cylinder collisions with walls.
#         '''
#         center_bound = self.grid_size - self.cyl_radii
#
#         left_collision = (self.config[:, 0] - self.cyl_radii < -self.grid_size)
#         right_collision = (self.config[:, 0] + self.cyl_radii > self.grid_size)
#
#         self.config[left_collision, 0] = -center_bound
#         self.config[right_collision, 0] = center_bound
#         self.vel[left_collision, 0] *= -self.wall_bounce
#         self.vel[right_collision, 0] *= -self.wall_bounce
#
#         top_collision = (self.config[:, 1] + self.cyl_radii > self.grid_size)
#         bottom_collision = (self.config[:, 1] - self.cyl_radii < -self.grid_size)
#
#         self.config[top_collision, 1] = center_bound
#         self.config[bottom_collision, 1] = -center_bound
#         self.vel[top_collision, 1] *= -self.wall_bounce
#         self.vel[bottom_collision, 1] *= -self.wall_bounce
