# ===== Dynamics and Simulations ===== #
# October 2026
# Esraaj Sarkar Gupta
#
#
# Ported from my old python codes for the same
# simulations.

# We simulate a double pendulum here using the Lagrangian Formalism.

using DifferentialEquations

include("system.jl")
include("animations.jl")

using .System
using .Animate

const OUTPUT_DIR = joinpath(@__DIR__, "figures/")
mkpath(OUTPUT_DIR)

params = DoublePenduParams(m1=1.0, m2=1.0, l1=1.0, l2=1.0, g=9.81)

# Initial conditions: [θ1, θ2, ω1, ω2]
u0 = [π / 2, π / 2, 0.0, 0.0]
tspan = (0.0, 20.0)

prob = ODEProblem(double_pendulum, u0, tspan, params)
sol = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-9)

# ===== Figure & Animation Generation via Animate Module ===== #

# 1. State Trajectories Plot
plot_trajectories(sol, OUTPUT_DIR)

# 2. Phase Space Plot
plot_phase_space(sol, OUTPUT_DIR)

# 3. End-Effector Spatial Trajectory (Cartesian)
plot_cartesian_path(sol, params, OUTPUT_DIR)

# 4. Animation Snapshot Screenshot
render_snapshot(sol, params, OUTPUT_DIR)

# 5. Full GIF Animation
animate_pendulum(sol, params, OUTPUT_DIR, fps=30)

println("Simulation & Visualization complete. Files saved to: ", OUTPUT_DIR)
