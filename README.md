# DesignEnvironments.jl

## Example of optimal control

<p>This is an example of an agent interacting with the CylinderEnv</p>

<p align="center">
<img src="https://github.com/AmirkulovaLab/DesignEnvironments.jl/blob/main/images/anim.gif" width="600">
</p>

## Example usage

```
using DesignEnvironments
using ReinforcementLearning
ENV["GKSwstype"] = "nul"
## need to set this ENV variable to nul in order to prevent animation from hanging

## create CylinderEnv object
env = CylinderEnv(M=10, continuous=true, grid_size=10.0)

## create a random policy suited to env action space
policy = RandomPolicy(action_space(env))

## generate video
render(env, policy, path="thing.mp4")
```
