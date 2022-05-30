var documenterSearchIndex = {"docs":
[{"location":"#DesignEnvironments.jl-Documentation","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.jl Documentation","text":"","category":"section"},{"location":"","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.jl Documentation","text":"                       _        _            _                      _             _     \n     /\\               (_)      | |          | |                    | |           | |    \n    /  \\    _ __ ___   _  _ __ | | __ _   _ | |  ___ __   __ __ _  | |      __ _ | |__  \n   / /\\ \\  | '_ ` _ \\ | || '__|| |/ /| | | || | / _ \\\\ \\ / // _` | | |     / _` || '_ \\\n  / ____ \\ | | | | | || || |   |   < | |_| || || (_) |\\ V /| (_| | | |____| (_| || |_) |\n /_/    \\_\\|_| |_| |_||_||_|   |_|\\_\\ \\__,_||_| \\___/  \\_/  \\__,_| |______|\\__,_||_.__/","category":"page"},{"location":"","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.jl Documentation","text":"","category":"page"},{"location":"","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.jl Documentation","text":"<style>\n.row {\n  display: flex;\n}\n\n.column {\n  flex: 50%;\n  padding: 5px;\n}\n</style>\n\n<div class=\"row\">\n    <div class=\"column\">\n        <img src=\"assets/config.gif\">\n        <p>Figure 1: A Configuration design which manages several cylindrical scatterers on a plane. The Total Scattering Cross Section is displayed on the right.</p>\n    </div>\n    <div class=\"column\">\n        <img src=\"assets/core_config.gif\">\n        <p>Figure 2: A CoreConfiguration design. This design has a static core which acts as an object which produces scattering. The mobile cylinders are controlled by the agent with the objective of shielding the core. The scattering produced by the entire configuration is shown in blue while the core scattering is shown in orange.</p>\n    </div>\n</div>","category":"page"},{"location":"#DesignEnvironment","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironment","text":"","category":"section"},{"location":"","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.jl Documentation","text":"<p>\nAn inverse design problem requires two components: a design, and an objective. In this package we provide an interface struct which holds these two components. The design encapsulates a physical system which can be modified through design actions. While the objective is a metric which is calculated on the design. This metric can be minimized or maximized according to the reward function.\n</p>","category":"page"},{"location":"","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.jl Documentation","text":"DesignEnvironment","category":"page"},{"location":"#DesignEnvironments.DesignEnvironment","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.DesignEnvironment","text":"Base environment for episodic design optimization.\n\nExample\n\nenv = DesignEnvironment(\n    design = Configuration(M = M, plane = Square(10.0), vel_decay=0.9),\n    objective = TSCS(k0amax = 1.0, k0amin = 0.3, nfreq = 15)\n    )\n\nParameters\n\ndesign::AbstractDesign\nobjective::AbstractObjective\nis_continuous::Bool: specifies if design actions are continuous or discrete\nepisode_length::Int: number of steps (actions) in a design episode\npenalty_weight::Float64: scalar which determines penalty for actions which result in invalid states\ntimestep::Int: current step number in the episode\npenalty::Float64: holds the current penalty\n\n\n\n\n\n","category":"type"},{"location":"#Designs","page":"DesignEnvironments.jl Documentation","title":"Designs","text":"","category":"section"},{"location":"","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.jl Documentation","text":"Configuration\nConfiguration(;::Int, ::DE.AbstractPlane, ::Float64, ::Float64, ::Float64)\nConfiguration(::Matrix, ::Matrix, ::Vector, ::DE.AbstractPlane, ::Float64, ::Float64, ::Float64)\nConfiguration(::Any)\nCoreConfiguration","category":"page"},{"location":"#DesignEnvironments.Configuration","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.Configuration","text":"Simulates a configuration of cylindrical scatterers on a two dimentional plane.\n\nParameters\n\nM::Int: number of cylindrical scatterers\nplane::AbstractPlane: 2 dimentional surface on which scatterers exist\nmax_vel::Float64: maximum allowable movement speed of cylinders\nvel_decay::Float64: percentage of velocity which remains after each step (0, 1.0)\nmin_distance::Float64: minimum allowable distance between cylinders\npos::Matrix: (x, y) coordinates of each cylinder\nvel::Matrix: (dx, dy) velocity of each cylinder\nradii::Vector: radius of each cylinder\n\n\n\n\n\n","category":"type"},{"location":"#DesignEnvironments.Configuration-Tuple{Matrix, Matrix, Vector, DesignEnvironments.AbstractPlane, Float64, Float64, Float64}","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.Configuration","text":"Constructor for defining a Configuration from a preexisting set of coordinates.\n\nExample\n\npos = zeros(M, 2)\nvel = zeros(M, 2)\nradii = ones(M)\n\nconfig = Configuration(\n    pos, vel, radii, Square(15.0), \n    0.1, 0.7, 0.1)\n\n\n\n\n\n","category":"method"},{"location":"#DesignEnvironments.Configuration-Tuple{Any}","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.Configuration","text":"Function which handles the application of an action to the design. The action can be either a matrix of velocity adjustments or an integer which represents a discrete action.\n\nThe function returns a penalty from the action which is to be minimized in the reward function.\n\nArguments\n\naction::Union{Vector, Int}\n\nReturns\n\n::Int\n\nExample\n\naction = randn(config.M * 2)\npenalty = config(action)\n\n\n\n\n\n","category":"method"},{"location":"#DesignEnvironments.CoreConfiguration","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.CoreConfiguration","text":"Applies an action to the Configuration of cylinders surrounding the core. Returns a penalty associated with the inter cylinder collisions and wall collisions.\n\n\n\n\n\n","category":"type"},{"location":"#Objectives","page":"DesignEnvironments.jl Documentation","title":"Objectives","text":"","category":"section"},{"location":"","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.jl Documentation","text":"TSCS\nTSCS(::Matrix)\nTSCS(::Configuration)\nTSCS(::CoreConfiguration)","category":"page"},{"location":"#DesignEnvironments.TSCS","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.TSCS","text":"Calculates the Total Scattering Cross Section of a planar Configuration of cylindrical scatterers. TSCS is calculated independently for each wavenumber (nfreq). Therefore this calculation benifits from multithreading. Start julia with threads <= nfreq to take advantage of this.\n\nExample\n\nobjective = TSCS(k0amax = 1.0, k0amin = 0.3, nfreq = 30)\n\nParameters\n\nk0amax::Float64\nk0amin::Float64\nnfreq::Int\nrho::Float64\nc0::Float64\naa::Float64\na::Float64\nha::Float64\nQ_RMS::Float64\nqV::Vector\nQ::Vector\n\n\n\n\n\n","category":"type"},{"location":"#DesignEnvironments.TSCS-Tuple{Configuration}","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.TSCS","text":"Defines function call for TSCS on a Configuration\n\nExample\n\nobjective(config)\n\n\n\n\n\n","category":"method"},{"location":"#DesignEnvironments.TSCS-Tuple{CoreConfiguration}","page":"DesignEnvironments.jl Documentation","title":"DesignEnvironments.TSCS","text":"Defines function call for TSCS on a CoreConfiguration\n\n\n\n\n\n","category":"method"}]
}
