# ===== Engineering Mathematics -- Homework 3, Question 8 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 15th February, 2026
#
# Repeated in Dynamics and Simulations
# Date: 20th August, 2026

mkpath("Problem_9/figures")

using DifferentialEquations
using LinearAlgebra
using Plots

# --- Define Parameters --- #
m   = 1
θ   = pi/4
v₀  = 10
g   = 1
c   = 1

k = c/m

t_final = 20

t = range(.0, t_final; length=1000)

# --- Initial Conditions --- #
x₀, y₀ = .0, .0
vx₀ = v₀ * cos(θ)
vy₀ = v₀ * sin(θ)

# ===== Part A ===== #
# I have solved for the system analytically, I may not be correct,
# but I have tried.

function projectile_analytic(t, vx₀, vy₀, k, g)
    x = (vx₀ / k) * (1 - exp(-k * t))
    y = (vy₀ / k) * (1 - exp(-k * t)) - (g / k) * (t - (1 / k) * (1 - exp(-k * t)))
    return x, y
end

results = projectile_analytic.(t, vx₀, vy₀, k, g)
x_analytic = [r[1] for r in results]
y_analytic = [r[2] for r in results]

idx_ground = findfirst(val -> val < 0, y_analytic[2:end]) + 1
stop_at = isnothing(idx_ground) ? length(y_analytic) : idx_ground

x_analytic_sliced = x_analytic[1:stop_at]
y_analytic_sliced = y_analytic[1:stop_at]

trajectory = plot(x_analytic_sliced, y_analytic_sliced, aspect_ratio=:equal,  xlabel="x", ylabel="y", label="analytic solution")

# ===== Part B ===== #

# State u = [x, y, dx, dy]
function projectile!(du, u, p, t)
    # Unpack parameter vector
    m, c, g = p

    k = c/m # For convinience

    # Unpack state vector
    x, y, dx, dy = u

    ddy = (- g) + (-k * dy)
    ddx = -k * dx

    du[1] = dx
    du[2] = dy
    du[3] = ddx
    du[4] = ddy
end

# -- Initial Conditions -- #

u₀  = [x₀, y₀, vx₀, vy₀]
tspan = (.0, t_final)

# Parameter Vector
p = (m, c, g)

# Callback conditon -- stop when y = 0 (hits the floor)

condition(u, t, integrator) = u[2]
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition, affect!)

problem = ODEProblem(projectile!, u₀, tspan, p)
sol = solve(problem, Tsit5(), callback = cb)

x_solutions = sol[1,:]
y_solutions = sol[2,:]
times = sol.t

plot!(trajectory, x_solutions, y_solutions,aspect_ratio=:equal, xlabel="x (m)", ylabel="y (m)", label="numerical solution")
savefig(trajectory, "Problem_9/figures/trajectory.png")

# ===== Part C ===== #

# Evaluate the analytical solution at (t=2)
state_at_2 = [projectile_analytic(2, vx₀, vy₀, k, g)...]
errors = []

tolerances = 10 .^ range(-3, -16, length=20)

for tol in tolerances
    temp_sol = solve(problem, Tsit5(), callback=cb, reltol=tol, abstol=tol)
    
    u_num = temp_sol(2.0)
    pos_num = [u_num[1], u_num[2]]
    
    push!(errors, norm(pos_num - state_at_2))
end

tol_error = plot(tolerances, errors, 
     xscale=:log10, yscale=:log10, 
     marker=:circle, 
     xlabel="Tolerance", ylabel="Error at t=2",
     title="Error vs Tolerance", flipx=true)

savefig(tol_error, "Problem_9/figures/tolerance_error.png")


