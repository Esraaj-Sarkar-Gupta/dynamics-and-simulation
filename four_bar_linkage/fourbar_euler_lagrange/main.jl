# Run: julia --project=. main.jl
include("symbolics.jl")
include("system.jl")
include("animations.jl")

using Printf
using LinearAlgebra
using .FourBarSystem
using .FourBarAnimations

function main()
    # Grashof crank-rocker geometry: theta1 can rotate through 2pi.
    params = FourBarParameters(
        l0=1.5, l1=0.6, l2=1.4, l3=1.0,
        m1=1.0, m2=1.2, m3=0.9, g=9.81
    )

    # Set torque=t -> some_expression to actuate link OA at O.
    # A zero torque gives freely evolving motion under gravity.
    model = build_system(params; torque=t -> 0.0)

    # theta2 and omega2 are initialized from the single length constraint.
    # branch=+1 or -1 selects the assembly mode for the INITIAL position.
    # Thereafter the system evolves both theta1 and theta2 dynamically.
    fps = 30
    sol = simulate(model;
                   theta1=pi/3, omega1=0.0, branch=+1,
                   tspan=(0.0, 8.0), fps=fps)

    paths = save_visualizations(model, sol; outdir="output", fps=fps)

    max_phi = maximum(abs(closure_residual(model, u[1:2])) for u in sol.u)
    max_length_error = maximum(
        abs(norm(body_positions(params, u[1:2]).B -
                 body_positions(params, u[1:2]).D) - params.l3)
        for u in sol.u
    )
    initial_energy = total_energy(model, sol.u[1])
    energy_drift = maximum(abs(total_energy(model, u) - initial_energy)
                           for u in sol.u)

    @printf("Integration retcode:          %s\n", string(sol.retcode))
    @printf("Saved frames:                 %d\n", length(sol.t))
    @printf("Max |Phi| (squared length):   %.3e m^2\n", max_phi)
    @printf("Max |DB length - l3|:        %.3e m\n", max_length_error)
    @printf("Max mechanical energy drift:  %.3e J (zero input torque)\n", energy_drift)
    println("Trajectories: ", paths.trajectories_path)
    println("Snapshot:     ", paths.snapshot_path)
    println("Animation:    ", paths.animation_path)
end

if true
    main()
    print("Done!")
end
