export DesignEnv

abstract type DesignEnv <: AbstractEnv end

## override ReinforcementLearning functions
RLBase.action_space(env::DesignEnv) = env.action_space
RLBase.state_space(env::DesignEnv) = env.state_space

RLBase.reward(env::DesignEnv) = DE.reward(env)
RLBase.is_terminated(env::DesignEnv) = DE.is_terminated(env)
RLBase.state(env::DesignEnv) = DE.state(env)
RLBase.reset!(env::DesignEnv) = DE.reset!(env)
(env::DesignEnv)(action) = env(action)

function render(
        env::DesignEnv, policy::AbstractPolicy;
        max_tscs::Float64=5.0, path::String="anim.mp4",
        fps::Int=30
        )

    reset!(env)

    dif = (env.k0amax - env.k0amin) / env.nfreq
    freqv = range(env.k0amin, env.k0amax, length=env.nfreq) |> collect

    a = Animation()
    prog = ProgressUnknown("Working hard:", spinner=true)

    while !is_terminated(env)
        ProgressMeter.next!(prog)

        ## this plot is the image which displays the current state of the environment
        plot_1 = img(env)

        ## this plot shows the scattering pattern (TSCS) produced by the current configuration
        plot_2 = plot(
            freqv, env.Q,
            xlabel="ka", ylabel="TSCS",
            xlim=(freqv[1], freqv[end]),
            ylim=(0, max_tscs))

        ## create a side by side plot containing two subplots
        big_plot = plot(
            plot_1,
            plot_2,
            layout=@layout([a{0.6w} b]),
            )

        ## ust this plot as a frame in the animation
        frame(a, big_plot)

        ## apply the policy's action to the environment
        env(policy(env))
    end
    ProgressMeter.finish!(prog)

    ## convert collections of images into gif
    gif(a, path, fps=fps)
    closeall()
end
