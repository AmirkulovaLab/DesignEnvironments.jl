# DesignEnvironments.jl

<h1>
_____          .__        __         .__                      .____          ___.    
/  _  \   _____ |__|______|  | ____ __|  |   _______  _______  |    |   _____ \_ |__  
/  /_\  \ /     \|  \_  __ \  |/ /  |  \  |  /  _ \  \/ /\__  \ |    |   \__  \ | __ \
/    |    \  Y Y  \  ||  | \/    <|  |  /  |_(  <_> )   /  / __ \|    |___ / __ \| \_\ \
\____|__  /__|_|  /__||__|  |__|_ \____/|____/\____/ \_/  (____  /_______ (____  /___  /
     \/      \/               \/                            \/        \/    \/    \/
</h1>

## Example of optimal control
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
