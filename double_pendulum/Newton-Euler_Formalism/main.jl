# ===== Double Pendulum Dynamics Pipeline ===== #
# Esraaj Sarkar Gupta
# November 2026
#
# Main Orchestrator Script -- Newton--Euler formulation

include("symbolics.jl")
include("system.jl")
include("animations.jl")

using .SymbolicSystem
using .System
using .Animations

function main()
    # 1. Specify output directory for all generated figures and animations.
    output_dir = joinpath(@__DIR__, "figures")
    mkpath(output_dir)
    println("Output Directory: ", output_dir)

    # 2. Symbolically derive the Newton--Euler equations about both joints.
    println("\n[1/6] Deriving the joint-moment equations symbolically...")
    print_symbolic_derivation()

    # 3. Define system parameters.
    params = DoublePendulumParams(
        m1 = 1.0, # kg
        m2 = 1.0, # kg
        l1 = 1.0, # m
        l2 = 1.0, # m
        g  = 9.81 # m/s^2
    )

    # 4. Define initial conditions: angles from the downward vertical.
    θ1_0, θ2_0 = π/2, π/2
    ω1_0, ω2_0 = 0.0, 0.0
    u0 = initial_conditions(θ1_0, θ2_0; ω1=ω1_0, ω2=ω2_0)
    tspan = (0.0, 10.0)

    # 5. Solve A(θ) * α = b(θ, ω), integrated as an ODE with Tsit5.
    println("\n[2/6] Solving Newton--Euler ODE via matrix formulation A*α = b...")
    sol = solve_double_pendulum(u0, tspan, params; saveat=0.02)
    println("      Solved $(length(sol.t)) saved time steps.")

    # 6. Generate all plots and animations through animations.jl.
    println("\n[3/6] Generating trajectory plot for both masses...")
    create_trajectory_plot(sol, params, output_dir; filename="trajectories.png")

    println("\n[4/6] Saving snapshot for the report...")
    save_snapshot(sol, params, output_dir; frame_idx=150,
                  filename="report_snapshot.png")

    println("\n[5/6] Checking numerical energy conservation...")
    create_energy_plot(sol, params, output_dir;
                       filename="energy_conservation.png")

    println("\n[6/6] Generating animated GIF...")
    create_animation(sol, params, output_dir;
                     filename="double_pendulum.gif", fps=30)

    println("\nAll outputs saved successfully in: $output_dir")
end

main()
