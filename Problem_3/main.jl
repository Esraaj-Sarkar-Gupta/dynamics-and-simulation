using DifferentialEquations
using Plots
using LinearAlgebra # Vector norms

mkpath("figures")

# Load dynamics and custom solver modules
include("dynamics.jl")
include("mySolvers.jl")
using .MySolvers

# System Parameters
u0 = [1.0, 1.0, 1.0]
p = [10.0, 28.0, 8/3]
tspan = (0.0, 20.0)

# ===========================
# Part (a): Direct Euler Loop
# ===========================
dt_direct = 0.001
N = round(Int, (tspan[2] - tspan[1]) / dt_direct)

x_dir = zeros(N + 1); y_dir = zeros(N + 1); z_dir = zeros(N + 1); t_dir = zeros(N + 1)

x_dir[1], y_dir[1], z_dir[1] = u0
t_dir[1] = tspan[1]

for i in 1:N
    dx = p[1] * (y_dir[i] - x_dir[i])
    dy = x_dir[i] * (p[2] - z_dir[i]) - y_dir[i]
    dz = x_dir[i] * y_dir[i] - p[3] * z_dir[i]

    x_dir[i + 1] = x_dir[i] + dt_direct * dx
    y_dir[i + 1] = y_dir[i] + dt_direct * dy
    z_dir[i + 1] = z_dir[i] + dt_direct * dz
    t_dir[i + 1] = t_dir[i] + dt_direct
end

plt_a = plot(x_dir, y_dir, z_dir, 
    xlabel="x", ylabel="y", zlabel="z", 
    title="Part (a): Direct Euler Loop", legend=false, lw=0.7
)
savefig(plt_a, "figures/lorenz_direct.png")


# ===================================================
# Parts (b) & (c): Custom Euler vs. ODE45 Comparison
# ===================================================
dt = 0.001

# Solve with custom Euler solver
sol_euler = MySolvers.Euler(Lorentz_system!, u0, tspan, p, dt)

# Solve with ODE45
prob = ODEProblem(Lorentz_system!, u0, tspan, p)
sol_ode45 = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# Unpack Euler state vectors
t_e = sol_euler.t
x_e = [u[1] for u in sol_euler.u]
y_e = [u[2] for u in sol_euler.u]
z_e = [u[3] for u in sol_euler.u]

# Interpolate ODE45 to match Euler's time grid
u_ode45_interp = sol_ode45(t_e)
x_o = [u[1] for u in u_ode45_interp]

# Plot 1: Overlay 3D Trajectories
plt_3d_comp = plot(x_e, y_e, z_e, label="Custom Euler (dt=0.001)", color=:red, lw=0.7)
plot!(plt_3d_comp, sol_ode45, idxs=(1,2,3), label="ODE45 (DP5)", color=:blue, lw=0.7,
      title="Lorenz: Custom Euler vs ODE45 Trajectory", xlabel="x", ylabel="y", zlabel="z")
savefig(plt_3d_comp, "figures/comparison_euler_vs_ode45_3d.png")

# Plot 2: Time Series Comparison x(t) (Highlighting Chaotic Divergence)
plt_time_comp = plot(t_e, x_e, label="Custom Euler", color=:red, lw=1.2)
plot!(plt_time_comp, t_e, x_o, label="ODE45", color=:blue, ls=:dash, lw=1.2,
      title="x(t) Trajectory Divergence (Butterfly Effect)", xlabel="Time (s)", ylabel="x(t)")
savefig(plt_time_comp, "figures/comparison_euler_vs_ode45_time.png")


# =================================================================
# Part (d): Euler Step-Size Convergence & Self-Accuracy Estimation
# =================================================================
# Estimating accuracy without an exact solution:

# We compute a high-precision reference solution using an 8th-order solver (Vern8)
sol_ref = solve(prob, Vern8(), reltol=1e-13, abstol=1e-13)

# Evaluate at T_eval = 5.0
T_eval = 5.0
u_ref_T = sol_ref(T_eval)

dt_range = [0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001]
errors_euler = Float64[]

for dt_val in dt_range
    sol_temp = MySolvers.Euler(Lorentz_system!, u0, (0.0, T_eval), p, dt_val)
    u_final_euler = sol_temp.u[end]
    push!(errors_euler, norm(u_final_euler - u_ref_T))
end

# Convergence Plot (Log-Log)
plt_conv = plot(dt_range, errors_euler, xscale=:log10, yscale=:log10, 
                marker=:circle, lw=2, label="Euler Error at T=5.0",
                xlabel="Step Size (dt)", ylabel="Absolute Error ||u_euler - u_ref||",
                title="Euler Method Convergence (Part d)")

# Overlay theoretical O(dt^1) reference slope
plot!(plt_conv, dt_range, dt_range .* (errors_euler[end] / dt_range[end]), 
      label="Theoretical O(dt¹)", ls=:dash, color=:black, lw=1.5)
savefig(plt_conv, "figures/euler_convergence.png")


# =====================================================
# Part (e): ODE45 Tolerance vs. Actual Error Assessment
# =====================================================
target_tols = 10.0 .^ (-3:-1:-10)
actual_errors = Float64[]

# High-precision ground truth at T = 10.0
T_tol_eval = 10.0
u_true_T = sol_ref(T_tol_eval)

for tol in target_tols
    sol_tol = solve(prob, DP5(), reltol=tol, abstol=tol)
    u_final_tol = sol_tol(T_tol_eval)
    push!(actual_errors, norm(u_final_tol - u_true_T))
end

# Tolerance Assessment Plot
plt_tol = plot(target_tols, actual_errors, xscale=:log10, yscale=:log10,
               marker=:square, lw=2, label="Actual Global Error",
               xlabel="Target Tolerance (reltol = abstol)", 
               ylabel="Actual Global Error at T=10.0",
               title="ODE45 Accuracy Self-Estimation (Part e)")

plot!(plt_tol, target_tols, target_tols, label="Ideal Error = Tolerance", 
      ls=:dash, color=:black, lw=1.5)
savefig(plt_tol, "figures/ode45_accuracy_assessment.png")
println("We're done here!")