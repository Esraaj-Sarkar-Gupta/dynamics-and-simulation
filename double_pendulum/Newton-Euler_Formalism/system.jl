# ===== Double Pendulum: Newton--Euler System ===== #
# Numerical equations of motion derived by taking moments about both joints.

module System

using LinearAlgebra
using DifferentialEquations

export DoublePendulumParams, construct_system_matrix, construct_rhs_vector,
       double_pendulum_ode!, solve_double_pendulum, initial_conditions,
       angles_to_cartesian, total_energy

Base.@kwdef struct DoublePendulumParams
    m1::Float64 = 1.0  # kg
    m2::Float64 = 1.0  # kg
    l1::Float64 = 1.0  # m
    l2::Float64 = 1.0  # m
    g::Float64  = 9.81 # m/s^2
end

"""
Return the 2×2 coefficient matrix A(θ) from the Newton--Euler equations,
written as A(θ) * [θ̈1, θ̈2] = b(θ, ω).

State ordering: u = [θ1, θ2, ω1, ω2].
"""
function construct_system_matrix(u, p::DoublePendulumParams)
    θ1, θ2 = u[1:2]
    Δ = θ1 - θ2
    c = cos(Δ)

    A = zeros(Float64, 2, 2)
    A[1, 1] = (p.m1 + p.m2) * p.l1
    A[1, 2] = p.m2 * p.l2 * c
    A[2, 1] = p.l1 * c
    A[2, 2] = p.l2

    return A
end

"""
Return the right-hand side b(θ, ω) of the reduced Newton--Euler system.

The first equation comes from the moment about O minus the moment about A;
the second comes from the moment about the accelerating joint A.
"""
function construct_rhs_vector(u, p::DoublePendulumParams)
    θ1, θ2, ω1, ω2 = u
    Δ = θ1 - θ2
    s = sin(Δ)

    b = zeros(Float64, 2)
    b[1] = -p.m2 * p.l2 * ω2^2 * s - (p.m1 + p.m2) * p.g * sin(θ1)
    b[2] =  p.l1 * ω1^2 * s - p.g * sin(θ2)

    return b
end

"""First-order right-hand side for DifferentialEquations.jl."""
function double_pendulum_ode!(du, u, p::DoublePendulumParams, t)
    du[1] = u[3]
    du[2] = u[4]

    A = construct_system_matrix(u, p)
    b = construct_rhs_vector(u, p)
    α = A \ b  # α = [θ̈1, θ̈2]

    du[3] = α[1]
    du[4] = α[2]
    return nothing
end

"""Construct angular initial conditions; all angles are from downward vertical."""
function initial_conditions(θ1, θ2; ω1=0.0, ω2=0.0)
    return Float64[θ1, θ2, ω1, ω2]
end

"""Convert joint angles into Cartesian positions of both point masses."""
function angles_to_cartesian(θ1, θ2, p::DoublePendulumParams)
    x1 = p.l1 * sin(θ1)
    y1 = -p.l1 * cos(θ1)
    x2 = x1 + p.l2 * sin(θ2)
    y2 = y1 - p.l2 * cos(θ2)
    return x1, y1, x2, y2
end

"""Total mechanical energy T + V (potential zero at y = 0)."""
function total_energy(u, p::DoublePendulumParams)
    θ1, θ2, ω1, ω2 = u
    Δ = θ1 - θ2

    T = 0.5 * (p.m1 + p.m2) * p.l1^2 * ω1^2 +
        0.5 * p.m2 * p.l2^2 * ω2^2 +
        p.m2 * p.l1 * p.l2 * ω1 * ω2 * cos(Δ)

    V = -(p.m1 + p.m2) * p.g * p.l1 * cos(θ1) -
        p.m2 * p.g * p.l2 * cos(θ2)

    return T + V
end

"""Integrate the Newton--Euler equations using Tsit5."""
function solve_double_pendulum(u0, tspan, p::DoublePendulumParams; saveat=0.02)
    (p.m1 > 0 && p.m2 > 0 && p.l1 > 0 && p.l2 > 0) ||
        throw(ArgumentError("Masses and rod lengths must be positive."))

    prob = ODEProblem(double_pendulum_ode!, u0, tspan, p)
    return solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat=saveat)
end

end # module System
