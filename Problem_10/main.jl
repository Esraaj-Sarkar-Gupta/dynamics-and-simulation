# ===== Dynamics and Simulations ===== #
# Author: Esraaj Sarkar Gupta
# Date: September, 2026
#

using LinearAlgebra
using DifferentialEquations
using Plots

# ---- Define the System ---- #

function system!(du, u, p, t)
    # Unpack Parameters
    m, k, L0, g = p

    # Views avoid memory allocation during integration
    @views r  = u[1:3]
    @views v  = u[4:6]
    @views dr = du[1:3]
    @views dv = du[4:6]

    # Kinematics
    dr .= v

    r_norm = norm(r)
    y_hat = [0.0, 0.0, 1.0]

    if r_norm > 0
        T = -k * (r_norm - L0) * (r / r_norm)
    else
        T = zero(r)
    end

    # Dynamics: a = T/m - g y_hat
    dv .= (T ./ m) .- (g .* y_hat)
end

# ---- Integrator Setup ---- #

# Parameters: m (kg), k (N/m), L0 (m), g (m/s^2)
p = (1.0, 5.0, 0.0, 9.81)

# Initial State: [x, y, z, vx, vy, vz]
# Mass pulled to (1, 1, 0.5) with initial velocity
u0 = [
    0.0, 0.0, 0.0,
    0.05, 0.05, 0.0
]
tspan = (0.0, 10.0)

# Solve ODE
prob = ODEProblem(system!, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# Extract spatial coordinates
x = sol[1, :]
y = sol[2, :]
z = sol[3, :]

# ---- 3D Trajectory Plot ---- #

plt = plot3d(
    x, y, z,
    label = "Trajectory",
    xlabel = "X (m)",
    ylabel = "Y (m)",
    zlabel = "Z (m)",
    title = "3D Spring-Mass Dynamics",
    lw = 1.5,
    color = :navy
)
scatter3d!([u0[1]], [u0[2]], [u0[3]], color=:green, ms=4, label="Start")
scatter3d!([x[end]], [y[end]], [z[end]], color=:red, ms=4, label="End")

savefig(plt, "Problem_10/figures/trajectory_3d.png")

# ---- Animation Generation ---- #

# Frame subsampling for smooth animation rendering
fps = 30
t_samples = range(tspan[1], tspan[2], length=Int(fps * tspan[2]))

anim = @animate for t in t_samples
    current_u = sol(t)
    rx, ry, rz = current_u[1], current_u[2], current_u[3]
    
    # Trace history up to time t
    past_sol = sol(range(tspan[1], t, length=200))
    px, py, pz = past_sol[1, :], past_sol[2, :], past_sol[3, :]

    # Plot origin anchor (spring root), spring line, and point mass
    plot3d(
        [0, rx], [0, ry], [0, rz],
        line = (:crimson, 2, :dash),
        label = "Spring",
        xlims = (minimum(x) - 0.2, maximum(x) + 0.2),
        ylims = (minimum(y) - 0.2, maximum(y) + 0.2),
        zlims = (minimum(z) - 0.2, maximum(z) + 0.2),
        xlabel = "X", ylabel = "Y", zlabel = "Z",
        title = "3D Spring-Mass System (t = $(round(t, digits=2))s)",
        legend = :topright
    )
    plot3d!(px, py, pz, color=:navy, alpha=0.5, lw=1.5, label="Path")
    scatter3d!([rx], [ry], [rz], color=:darkorange, ms=6, label="Mass (m)")
    scatter3d!([0], [0], [0], color=:black, ms=4, label="Anchor")
end

gif(anim, "Problem_10/figures/spring_mass_3d.gif", fps=fps)