# Four-Bar Linkage
# Esraaj Sarkar Gupta
include("symbolics.jl")
include("system.jl")
include("animations.jl")

using LinearAlgebra
using Printf
using .FourBarDAESystem
using .FourBarDAEAnimations

function main()
    params = FourBarParameters(
        l0=1.5, l1=0.6, l2=1.4, l3=1.0,
        m1=1.0, m2=1.2, m3=0.9, g=9.81,
    )
    model = build_system(params; torque=t -> 0.0)
    fps = 30
    sol = simulate(model; theta1=pi/3, omega1=0.0, branch=1,
                   tspan=(0.0,8.0), fps)
    outputs = save_visualizations(model, sol; outdir=joinpath(@__DIR__,"output"), fps)

    max_psi = maximum(norm(constraint_residual(model,u[1:4]),Inf) for u in sol.u)
    max_length_error = maximum(begin
        s=configuration(model,u)
        maximum(abs.([norm(s.rA)-params.l1,
                      norm(s.rB-s.rA)-params.l2,
                      norm(s.rB-s.rD)-params.l3]))
    end for u in sol.u)
    initial_energy=total_energy(model,sol.u[1])
    energy_drift=maximum(abs(total_energy(model,u)-initial_energy) for u in sol.u)

    @printf("Integration retcode:             %s\n",string(sol.retcode))
    @printf("Saved frames:                    %d\n",length(sol.t))
    @printf("Max |Psi|:                       %.3e m²\n",max_psi)
    @printf("Max link-length error:           %.3e m\n",max_length_error)
    @printf("Max mechanical energy drift:     %.3e J (zero input torque)\n",energy_drift)
    println("Trajectories: ",outputs.trajectories_path)
    println("Snapshot:     ",outputs.snapshot_path)
    println("Animation:    ",outputs.animation_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
