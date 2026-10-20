# ===== Dynamics and Simulations ===== #
# October 2026
# Esraaj Sarkar Gupta


# We simulate a double pendulum here using DAEs.

using DifferentialEquations
using LinearAlgebra

include("system.jl")
include("animations.jl")

using .System
using .Animations

function main()
    # Configuration
    save_dir = "double_pendulum/DAEs/figures/"
    
    if !isdir(save_dir) 
        mkpath(save_dir)
        println("Created directory: $save_dir")
    else
        println("Output directory exists: $save_dir")
    end

    # Instantiate parameters
    p = DoublePenduParams(m1=1.0, m2=1.0, l1=1.0, l2=1.0, g=9.81)

    # Initial conditions (must satisfy geometric constraints)
    θ1 = π / 2  # Horizontal right
    θ2 = π / 2  
    
    x1_0 = p.l1 * sin(θ1)
    y1_0 = -p.l1 * cos(θ1)
    
    x2_0 = x1_0 + p.l2 * sin(θ2)
    y2_0 = y1_0 - p.l2 * cos(θ2)

    # Initial velocities
    vx1_0, vy1_0 = 0.0, 0.0
    vx2_0, vy2_0 = 0.0, 0.0

    u0 = [x1_0, y1_0, x2_0, y2_0, vx1_0, vy1_0, vx2_0, vy2_0]
    tspan = (0.0, 10.0)

    println("Initializing solver...")
    prob = ODEProblem(double_pendulum!, u0, tspan, p)
    
    # Tsit5 is ideal for explicit first-order transformations 
    sol = solve(prob, Tsit5(), saveat=0.05, reltol=1e-8, abstol=1e-8)
    
    println("Integration complete. Generating outputs...")
    plot_trajectory(sol, save_dir)
    animate_pendulum(sol, p, save_dir, fps=30)
    
    println("All tasks finished successfully.")
end

main()