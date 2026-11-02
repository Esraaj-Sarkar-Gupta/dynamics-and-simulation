# ===== Triple Pendulum Dynamics Pipeline ===== #
# Newton--Euler formulation, taking moments about three joints.
# Main orchestrator script.

include("symbolics.jl")
include("system.jl")
include("animations.jl")

using .SymbolicSystem
using .System
using .Animations

function main()
    # 1. All figures go into a child directory.
    output_dir = joinpath(@__DIR__, "figures")
    mkpath(output_dir)
    println("Output Directory: ", output_dir)

    # 2. Derive equations rather than introducing them as unexplained formulas.
    #    The numeric system.jl contains the compact, equivalent coefficients.
    println("\n[1/7] Deriving and validating Newton--Euler equations symbolically...")
    print_symbolic_derivation()

    # 3. System parameters (same defaults as the Lagrangian/DAE examples).
    params = TriplePendulumParams(
        m1 = 1.0, m2 = 1.0, m3 = 1.0,  # kg
        l1 = 1.0, l2 = 1.0, l3 = 1.0,  # m
        g = 9.81                        # m/s^2
    )

    # 4. Angular initial conditions, measured from downward vertical.
    θ1_0, θ2_0, θ3_0 = π/2, π/2, π/2
    ω1_0, ω2_0, ω3_0 = 0.0, 0.0, 0.0
    u0 = initial_conditions(θ1_0, θ2_0, θ3_0;
                            ω1=ω1_0, ω2=ω2_0, ω3=ω3_0)
    tspan = (0.0, 10.0)

    # 5. Numerically solve A(θ) * α = b(θ,ω) with Tsit5.
    println("\n[2/7] Solving triple pendulum ODE...")
    sol = solve_triple_pendulum(u0, tspan, params; saveat=0.02)
    println("      Saved $(length(sol.t)) time steps.")

    # 6. All visualization and animation functions reside in animations.jl.
    println("\n[3/7] Generating angle plot...")
    create_angle_plot(sol, params, output_dir; filename="angles_plot.png")

    println("\n[4/7] Generating trajectories...")
    create_trajectory_plot(sol, params, output_dir; filename="trajectories.png")

    println("\n[5/7] Saving report snapshot...")
    save_snapshot(sol, params, output_dir;
                  frame_idx=150, filename="report_snapshot.png")

    println("\n[6/7] Checking energy conservation...")
    create_energy_plot(sol, params, output_dir;
                       filename="energy_conservation.png")

    println("\n[7/7] Generating animation GIF...")
    create_animation(sol, params, output_dir;
                     filename="triple_pendulum.gif", fps=30, frame_step=2)

    println("\nAll outputs saved to: ", output_dir)
    return sol
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
