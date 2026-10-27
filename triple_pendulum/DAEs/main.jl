# ===== Triple Pendulum Dynamics Pipeline ===== #
# Esraaj Sarkar Gupta
# October 2026
#
# Main Orchestrator Script -- DAEs

include("system.jl")
include("animations.jl")

using .System
using .Animations

function main()
    # 1. Specify child output directory for all generated files
    output_dir = joinpath(@__DIR__, "figures/")
    println("Output Directory: ", output_dir)

    # 2. Define System Parameters
    params = TriplePendulumParams(
        m1 = 1.0,  # kg
        m2 = 1.0,  # kg
        m3 = 1.0,  # kg
        l1 = 1.0,  # m
        l2 = 1.0,  # m
        l3 = 1.0,  # m
        g  = 9.81  # m/s^2
    )

    # 3. Define Initial Conditions (θ1, θ2, θ3 in radians)
    θ1_0, θ2_0, θ3_0 = π/2, π/2, π/2
    u0 = angles_to_cartesian(θ1_0, θ2_0, θ3_0, params)

    tspan = (0.0, 10.0)

    # 4. Numerical Solution using system.jl
    println("\n[1/4] Solving numerical ODE system via matrix formulation A*x = b...")
    sol = solve_triple_pendulum(u0, tspan, params; saveat=0.02)
    println("      Solved $(length(sol.t)) time steps.")

    # 5. Generate Outputs using animations.jl
    println("\n[2/4] Generating trajectory plot for all 3 end effectors...")
    create_trajectory_plot(sol, params, output_dir; filename="trajectories.png")

    println("\n[3/4] Saving snapshot frame for report...")
    save_snapshot(sol, params, output_dir; frame_idx=150, filename="report_snapshot.png")

    println("\n[4/4] Generating pendulum animation GIF...")
    create_animation(sol, params, output_dir; filename="triple_pendulum.gif", fps=30)

    println("\nAll outputs saved successfully in: $output_dir")
end

main()