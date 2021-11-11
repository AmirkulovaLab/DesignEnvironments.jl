```
                       _        _            _                      _             _     
     /\               (_)      | |          | |                    | |           | |    
    /  \    _ __ ___   _  _ __ | | __ _   _ | |  ___ __   __ __ _  | |      __ _ | |__  
   / /\ \  | '_ ` _ \ | || '__|| |/ /| | | || | / _ \\ \ / // _` | | |     / _` || '_ \
  / ____ \ | | | | | || || |   |   < | |_| || || (_) |\ V /| (_| | | |____| (_| || |_) |
 /_/    \_\|_| |_| |_||_||_|   |_|\_\ \__,_||_| \___/  \_/  \__,_| |______|\__,_||_.__/
```
# DesignEnvironments.jl

## Demo
<p align="center">
<img src="https://github.com/AmirkulovaLab/DesignEnvironments.jl/blob/main/images/physics.gif" width="600">
</p>

<p>This is an example of an agent interacting with the CylinderEnv. We can see that the agent is able to supress the TSCS across a range of wavenumbers.</p>

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
