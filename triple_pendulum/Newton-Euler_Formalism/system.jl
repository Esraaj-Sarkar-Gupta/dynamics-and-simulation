# ===== Triple Pendulum: Newton--Euler System ===== #
# Numeric dynamics obtained by taking moments about all three joints.
# The joint-moment derivation and extraction of A and b are in symbolics.jl.

module System

using LinearAlgebra
using DifferentialEquations

export TriplePendulumParams, initial_conditions, angles_to_cartesian,
       construct_system_matrix, construct_rhs_vector,
       triple_pendulum_ode!, solve_triple_pendulum, total_energy

Base.@kwdef struct TriplePendulumParams
    m1::Float64 = 1.0  # kg
    m2::Float64 = 1.0  # kg
    m3::Float64 = 1.0  # kg
    l1::Float64 = 1.0  # m
    l2::Float64 = 1.0  # m
    l3::Float64 = 1.0  # m
    g::Float64  = 9.81 # m/s^2
end

"""State vector [θ1, θ2, θ3, ω1, ω2, ω3], angles from downward vertical."""
function initial_conditions(θ1, θ2, θ3; ω1=0.0, ω2=0.0, ω3=0.0)
    return Float64[θ1, θ2, θ3, ω1, ω2, ω3]
end

"""Cartesian coordinates (x1,y1,x2,y2,x3,y3) for a given set of angles."""
function angles_to_cartesian(θ1, θ2, θ3, p::TriplePendulumParams)
    x1 = p.l1 * sin(θ1)
    y1 = -p.l1 * cos(θ1)
    x2 = x1 + p.l2 * sin(θ2)
    y2 = y1 - p.l2 * cos(θ2)
    x3 = x2 + p.l3 * sin(θ3)
    y3 = y2 - p.l3 * cos(θ3)
    return x1, y1, x2, y2, x3, y3
end

"""
Mass matrix A(θ) in A(θ) * [θ̈1, θ̈2, θ̈3] = b(θ, ω).
The three rows are obtained by subtracting successive Newton--Euler
joint-moment equations (E1-E2, E2-E3, E3).
"""
function construct_system_matrix(u, p::TriplePendulumParams)
    θ1, θ2, θ3 = u[1:3]
    l1, l2, l3 = p.l1, p.l2, p.l3
    S1 = p.m1 + p.m2 + p.m3
    S2 = p.m2 + p.m3
    S3 = p.m3

    A = zeros(Float64, 3, 3)
    A[1, 1] = S1 * l1^2
    A[1, 2] = A[2, 1] = S2 * l1 * l2 * cos(θ1 - θ2)
    A[1, 3] = A[3, 1] = S3 * l1 * l3 * cos(θ1 - θ3)
    A[2, 2] = S2 * l2^2
    A[2, 3] = A[3, 2] = S3 * l2 * l3 * cos(θ2 - θ3)
    A[3, 3] = S3 * l3^2
    return A
end

"""
Forcing vector from the Newton--Euler moments, after eliminating joint forces:
  b_i = -g*S_i*l_i*sin(θ_i)
        - Σ_{j≠i} S_max(i,j)*l_i*l_j*sin(θ_i-θ_j)*ω_j^2.
"""
function construct_rhs_vector(u, p::TriplePendulumParams)
    θ1, θ2, θ3, ω1, ω2, ω3 = u
    l1, l2, l3 = p.l1, p.l2, p.l3
    S1 = p.m1 + p.m2 + p.m3
    S2 = p.m2 + p.m3
    S3 = p.m3

    Δ12 = θ1 - θ2
    Δ13 = θ1 - θ3
    Δ23 = θ2 - θ3

    b1 = -S1*p.g*l1*sin(θ1) - S2*l1*l2*ω2^2*sin(Δ12) - S3*l1*l3*ω3^2*sin(Δ13)
    b2 = -S2*p.g*l2*sin(θ2) + S2*l1*l2*ω1^2*sin(Δ12) - S3*l2*l3*ω3^2*sin(Δ23)
    b3 = -S3*p.g*l3*sin(θ3) + S3*l1*l3*ω1^2*sin(Δ13) + S3*l2*l3*ω2^2*sin(Δ23)

    return Float64[b1, b2, b3]
end

"""First-order ODE right-hand side: u̇ = [ω1,ω2,ω3,α1,α2,α3]."""
function triple_pendulum_ode!(du, u, p::TriplePendulumParams, t)
    du[1:3] .= u[4:6]
    A = construct_system_matrix(u, p)
    b = construct_rhs_vector(u, p)
    du[4:6] .= A \ b  # Do not explicitly compute inv(A).
    return nothing
end

"""Compute T+V directly from Cartesian bob velocities and heights."""
function total_energy(u, p::TriplePendulumParams)
    θ1, θ2, θ3, ω1, ω2, ω3 = u
    l = (p.l1, p.l2, p.l3)
    m = (p.m1, p.m2, p.m3)
    θ = (θ1, θ2, θ3)
    ω = (ω1, ω2, ω3)

    xdot, ydot, y = 0.0, 0.0, 0.0
    T, V = 0.0, 0.0
    for i in 1:3
        xdot += l[i] * cos(θ[i]) * ω[i]
        ydot += l[i] * sin(θ[i]) * ω[i]
        y -= l[i] * cos(θ[i])
        T += 0.5 * m[i] * (xdot^2 + ydot^2)
        V += m[i] * p.g * y
    end
    return T + V
end

"""Numerically integrate the Newton--Euler system with Tsit5."""
function solve_triple_pendulum(u0, tspan, p::TriplePendulumParams;
                              saveat=0.02, reltol=1e-8, abstol=1e-8)
    all(>(0), (p.m1, p.m2, p.m3, p.l1, p.l2, p.l3)) ||
        throw(ArgumentError("All masses and rod lengths must be positive."))
    length(u0) == 6 || throw(ArgumentError("Initial state must contain 6 entries."))

    prob = ODEProblem(triple_pendulum_ode!, Float64.(u0), tspan, p)
    return solve(prob, Tsit5(); reltol=reltol, abstol=abstol, saveat=saveat)
end

end # module System
