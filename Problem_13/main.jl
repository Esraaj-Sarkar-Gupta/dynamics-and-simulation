# ===== Engineering Mathematics -- Homework 3, Question 9 ===== #
# Author: Esraaj Sarkar Gupta
# Date: 15th February, 2026

using DifferentialEquations
using LinearAlgebra
using Plots

# ---- Fundamental System Parameters ---- #

g = 9.810

# ---- Fundamental System Parameters ---- #

struct fparameters
    m  :: Float64
    rₐ :: Vector{Float64}
    rᵦ :: Vector{Float64}

    Lₐ :: Float64 
    Lᵦ :: Float64 
    
    kₐ :: Float64
    kᵦ :: Float64

    cₐ :: Float64
    cᵦ :: Float64

    c_d :: Float64
end

# ---- Derived System Parameters ---- #

struct dparameters
    lₐ :: Float64
    lᵦ :: Float64

    a :: Vector{Float64}
    b :: Vector{Float64}

end

function get_dparameters(sys::fparameters, r_c::Vector{Float64})
    a_vec = sys.rₐ - r_c
    b_vec = sys.rᵦ - r_c
    
    # Current scalar lengths
    lₐ_val = norm(a_vec)
    lᵦ_val = norm(b_vec)
    
    return dparameters(lₐ_val, lᵦ_val, a_vec, b_vec)
end

# ---- Define System Evolution ---- #

function spring_damper_mass_system!(du, u, p, t)
    # Unpack state
    r_c = u[1:2]
    v_c = u[3:4]

    # Unpack Parameters
    m = p.m 
    r_A, r_B = p.rₐ, p.rᵦ

    k_A, k_B = p.kₐ, p.kᵦ
    L_A, L_B = p.Lₐ, p.Lᵦ

    c_A, c_B, c_d = p.cₐ, p.cᵦ, p.c_d

    # -- Geometry -- #
    vec_to_A = r_A - r_c
    vec_to_B = r_B - r_c

    dist_A = norm(vec_to_A)
    dist_B = norm(vec_to_B)

    # Unit Vectors
    hat_A = vec_to_A / dist_A
    hat_B = vec_to_B / dist_B

    # -- Force Calculations -- #

    # - Springs - #
    delta_A = dist_A - L_A
    delta_B = dist_B - L_B
    
    F_spring_A = k_A * delta_A * hat_A
    F_spring_B = k_B * delta_B * hat_B

    # - Dampers - #
    vel_proj_A = dot(v_c, -hat_A)
    vel_proj_B = dot(v_c, -hat_B)

    F_damp_A = -c_A * dot(v_c, hat_A) * hat_A 
    F_damp_B = -c_B * dot(v_c, hat_B) * hat_B
    
    # - Drag - #
    F_drag = -c_d * v_c

    # - Gravity - #
    F_gravity = [0.0, -m * g]

    # -- Newton's Second Law -- #

    F_total = F_spring_A + F_spring_B + F_damp_A + F_damp_B + F_drag + F_gravity
    accel = F_total / m

    # -- Update State -- #
    du[1:2] = v_c    # dr/dt = v
    du[3:4] = accel  # dv/dt = F/m

end
# ---- Solving the Unremarkable System ---- #

params = fparameters(
    5.0,          # Mass (kg)
    [2.0, 10.0],  # Anchor A
    [10.0, 10.0], # Anchor B
    5.0,          # Rest Length A
    5.0,          # Rest Length B
    10.0,         # k_A
    20.0,         # k_B
    10.0,          # c_A
    1.0,          # c_B
    0.5           # c_d (Air Drag)
)

# Initial Conditions
u0 = [5.0, 0.0, 0.0, 2.0] 

# Time Span
tspan = (0.0, 50.0)

prob = ODEProblem(spring_damper_mass_system!, u0, tspan, params)
sol = solve(prob, Tsit5(), saveat=0.1)

x = sol[1, :]
y = sol[2, :]

show_path = true

# -- Convinience -- #

rA = params.rₐ
rB = params.rᵦ

LA = params.Lₐ * rA / norm(rA)

# ---- Padding and XY Limits ---- #

pad = 1.0
xmin = min(minimum(x), rA[1], rB[1]) - pad
xmax = max(maximum(x), rA[1], rB[1]) + pad
ymin = min(minimum(y), rA[2], rB[2]) - pad
ymax = max(maximum(y), rA[2], rB[2]) + pad

animation = @animate for i in eachindex(sol.t)
    xc, yc = x[i], y[i]

    p = plot(
        xlim=(xmin, xmax), ylim=(ymin, ymax),
        aspect_ratio=:equal,
        xlabel="x (m)", ylabel="y (m)",
        title="Mass + Springs    t = $(round(sol.t[i], digits=2)) s",
        legend=false
    )

    # Show background Path 
    if show_path
        plot!(p, x[1:i], y[1:i], lw=1, alpha=0.35)
    end

    # Springs
    plot!(p, [rA[1], xc], [rA[2], yc], lw=3)
    plot!(p, [rB[1], xc], [rB[2], yc], lw=3)

    # Anchors
    scatter!(p,
        [rA[1], rB[1]],
        [rA[2], rB[2]],
        markersize=8,
        markerstrokewidth=0
        )

    # Mass point
    scatter!(p, [xc], [yc], markersize=10, markerstrokewidth=0)

    p
end

gif(animation, "mass_springs.gif", fps=30)
#display(sol)