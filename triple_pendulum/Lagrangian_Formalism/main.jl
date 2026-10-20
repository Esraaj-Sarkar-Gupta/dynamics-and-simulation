# ===== Triple Pendulum Dynamics Pipeline ===== #
# Esraaj Sarkar Gupta
# October 2026
#
# Main Orchestrator Script

using LinearAlgebra
using DifferentialEquations

# Target child directory for output files
const OUTPUT_DIR = joinpath(@__DIR__, "figures/")

# Load modular components
include("symbolics.jl")
include("animations.jl")

# Build compiled Mass matrix and Forcing functions
M_func, f_func = derive_triple_pendulum_equations()

# Define ODE In-Place Dynamics
function triple_pendulum_dynamics!(du, u, p, t)
    M = M_func(u, p)
    f = f_func(u, p)
    
    # Fast numerical linear solve M \ f for angular accelerations
    ddth = M \ f
    
    du[1:3] .= u[4:6]
    du[4:6] .= ddth
end

# System parameters: [m1, m2, m3, l1, l2, l3, g]
p = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 9.81]

# Initial conditions: θ = π/2, ω = 0
u0 = [π/2, π/2, π/2, 0.0, 0.0, 0.0]

# Time span: 10 seconds
tspan = (0.0, 10.0)

println("4. Integrating ODE over t ∈ [0, 10] using Tsit5...")
prob = ODEProblem(triple_pendulum_dynamics!, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# Render and export graphics to child directory
plot_and_animate(sol, p, OUTPUT_DIR; snapshot_time=2.5)

println("\nAll tasks completed successfully. Outputs directory: ", OUTPUT_DIR)