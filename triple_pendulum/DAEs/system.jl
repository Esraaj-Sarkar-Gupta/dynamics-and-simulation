# system.jl
# Defines Newton's Second Law of Motion constraint_derivatives

module System

using LinearAlgebra
using DifferentialEquations

export TriplePendulumParams, solve_triple_pendulum, angles_to_cartesian

struct TriplePendulumParams
    m1::Float64
    m2::Float64
    m3::Float64
    l1::Float64
    l2::Float64
    l3::Float64
    g::Float64
end

# Default parameter constructor
TriplePendulumParams(; m1=1.0, m2=1.0, m3=1.0, l1=1.0, l2=1.0, l3=1.0, g=9.81) = 
    TriplePendulumParams(m1, m2, m3, l1, l2, l3, g)

"""
Constructs the 9x9 matrix A from the constrained formulation.
"""
function construct_system_matrix(u, p::TriplePendulumParams)
    x1, y1, x2, y2, x3, y3 = u[1:6]
    m1, m2, m3 = p.m1, p.m2, p.m3

    Δx2 = x2 - x1
    Δy2 = y2 - y1
    Δx3 = x3 - x2
    Δy3 = y3 - y2

    A = zeros(Float64, 9, 9)

    # 1. Mass Matrix Block (Rows 1-6, Cols 1-6)
    A[1, 1] = m1
    A[2, 2] = m1
    A[3, 3] = m2
    A[4, 4] = m2
    A[5, 5] = m3
    A[6, 6] = m3

    # 2. Jacobian Transpose J^T Block (Rows 1-6, Cols 7-9)
    A[1, 7] = x1;        A[1, 8] = -Δx2
    A[2, 7] = y1;        A[2, 8] = -Δy2
    A[3, 8] = Δx2;       A[3, 9] = -Δx3
    A[4, 8] = Δy2;       A[4, 9] = -Δy3
    A[5, 9] = Δx3
    A[6, 9] = Δy3

    # 3. Jacobian J Block (Rows 7-9, Cols 1-6)
    A[7, 1] = x1;        A[7, 2] = y1
    A[8, 1] = -Δx2;      A[8, 2] = -Δy2;      A[8, 3] = Δx2;   A[8, 4] = Δy2
    A[9, 3] = -Δx3;      A[9, 4] = -Δy3;      A[9, 5] = Δx3;   A[9, 6] = Δy3

    return A
end

"""
Constructs the 9x1 right-hand side vector b.
"""
function construct_rhs_vector(u, p::TriplePendulumParams)
    dx1, dy1, dx2, dy2, dx3, dy3 = u[7:12]
    m1, m2, m3, g = p.m1, p.m2, p.m3, p.g

    Δdx2 = dx2 - dx1
    Δdy2 = dy2 - dy1
    Δdx3 = dx3 - dx2
    Δdy3 = dy3 - dy2

    b = zeros(Float64, 9)
    b[1] = 0.0
    b[2] = -m1 * g
    b[3] = 0.0
    b[4] = -m2 * g
    b[5] = 0.0
    b[6] = -m3 * g
    b[7] = -(dx1^2 + dy1^2)
    b[8] = -(Δdx2^2 + Δdy2^2)
    b[9] = -(Δdx3^2 + Δdy3^2)

    return b
end

"""
ODE function evaluating du/dt = [velocities; accelerations]
"""
function triple_pendulum_ode!(du, u, p::TriplePendulumParams, t)
    # Velocity state derivatives
    du[1:6] .= u[7:12]

    # Solve A * x = b for accelerations and multipliers
    A = construct_system_matrix(u, p)
    b = construct_rhs_vector(u, p)
    x_sol = A \ b  # x_sol = [ddx1, ddy1, ddx2, ddy2, ddx3, ddy3, λ1, λ2, λ3]

    # Acceleration state derivatives
    du[7:12] .= x_sol[1:6]
    return nothing
end

"""
Utility function converting initial pendulum angles (θ1, θ2, θ3) relative 
to vertical into initial Cartesian positions and zero velocities.
"""
function angles_to_cartesian(θ1, θ2, θ3, p::TriplePendulumParams)
    x1 = p.l1 * sin(θ1)
    y1 = -p.l1 * cos(θ1)

    x2 = x1 + p.l2 * sin(θ2)
    y2 = y1 - p.l2 * cos(θ2)

    x3 = x2 + p.l3 * sin(θ3)
    y3 = y2 - p.l3 * cos(θ3)

    return [x1, y1, x2, y2, x3, y3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
end

"""
Solves the ODE system numerically using Tsit5 solver.
"""
function solve_triple_pendulum(u0, tspan, p::TriplePendulumParams; saveat=0.02)
    prob = ODEProblem(triple_pendulum_ode!, u0, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat=saveat)
    return sol
end

end # module System