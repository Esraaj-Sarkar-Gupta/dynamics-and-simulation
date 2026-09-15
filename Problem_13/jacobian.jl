# ===== Dynamics and Simulations -- Question 13 ===== #
# Author: Esraaj Sarkar Gupta
# September 2026
#
# Adapted partially from Robotic Simulation and Motion

using ForwardDiff
using NLsolve
using DifferentialEquations
using LinearAlgebra
using Plots

# ---- Fundamental System Parameters ---- #

g = 9.810

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

# ---- Solving Static System Forces ---- #

function static_forces!(F, r_c, sys::fparameters)
    a_vec = sys.rₐ - r_c
    b_vec = sys.rᵦ - r_c

    lₐ = norm(a_vec)
    lᵦ = norm(b_vec)

    hat_a = a_vec / lₐ
    hat_b = b_vec / lᵦ

    delta_a = lₐ - sys.Lₐ
    delta_b = lᵦ - sys.Lᵦ

    gravity = [0.0, -sys.m * g]
    force_a = sys.kₐ * delta_a * hat_a
    force_b = sys.kᵦ * delta_b * hat_b

    F[1] = gravity[1] + force_a[1] + force_b[1]
    F[2] = gravity[2] + force_a[2] + force_b[2]
end

# ---- System Parameters ---- #

params = fparameters(
    5.0,          # Mass (kg)
    [2.0, 10.0],  # Anchor A
    [10.0, 10.0], # Anchor B
    5.0,          # Rest Length A
    5.0,          # Rest Length B
    10.0,         # k_A
    20.0,         # k_B
    10.0,         # c_A
    1.0,          # c_B
    0.5           # c_d (Air Drag)
)

# Initial Conditions & Time Span
u0 = [5.0, 0.0, 0.0, 2.0] 
tspan = (0.0, 50.0)

prob = ODEProblem(spring_damper_mass_system!, u0, tspan, params)
sol = solve(prob, Tsit5(), saveat=0.1)

x = sol[1, :]
y = sol[2, :]

# ---- Static Root Finding & Custom Clustering ---- #

num_samples = 2000
solutions_x = Float64[]
solutions_y = Float64[]

target_func!(F, pos) = static_forces!(F, pos, params)

for i in 1:num_samples
    x0 = -5.0 + 20.0 * rand()
    y0 = -10.0 + 25.0 * rand()
    initial_guess = [x0, y0]

    res = nlsolve(target_func!, initial_guess, autodiff=:forward)

    if converged(res)
        push!(solutions_x, res.zero[1])
        push!(solutions_y, res.zero[2])
    end
end

# Combine points and form integer seed classes
points = tuple.(solutions_x, solutions_y)
rounded_points = tuple.(round.(Int, solutions_x), round.(Int, solutions_y))
unique_classes = unique(rounded_points)

println("Identified $(length(unique_classes)) integer class seeds.")

clusters = Dict(c => Tuple{Float64, Float64}[] for c in unique_classes)

# Assign points to nearest integer class using Euclidean distance
for p in points
    best_class = unique_classes[1]
    min_dist = Inf

    for c in unique_classes
        dist = sqrt((p[1] - c[1])^2 + (p[2] - c[2])^2)
        if dist < min_dist
            min_dist = dist
            best_class = c
        end
    end
    push!(clusters[best_class], p)
end

# Compute cluster centroids and construct 4D state vectors u* = [cx, cy, 0, 0]
fixed_points = Vector{Vector{Float64}}()

for (class_seed, cluster_points) in clusters
    if !isempty(cluster_points)
        cx = sum(p[1] for p in cluster_points) / length(cluster_points)
        cy = sum(p[2] for p in cluster_points) / length(cluster_points)
        
        # Fixed point in 4D state-space: zero velocity at static equilibrium
        push!(fixed_points, [cx, cy, 0.0, 0.0])
        println("Region Centroid: (X: $(round(cx, digits=4)), Y: $(round(cy, digits=4)))")
    end
end

# ---- Linearization & Eigenvalue Stability Analysis ---- #

function rhs(u)
    du = Vector{eltype(u)}(undef, 4)
    spring_damper_mass_system!(du, u, params, 0.0)
    return du
end

println("\n" * "="^50)
println("SYSTEM STABILITY ANALYSIS")
println("="^50)

for (i, fp) in enumerate(fixed_points)
    # Compute 4x4 Jacobian Matrix at u* via Forward-Mode AD
    J = ForwardDiff.jacobian(rhs, fp)
    
    # Calculate Eigenvalues
    λ = eigvals(J)
    
    # Stability criterion: all Re(λ) < 0
    is_stable = all(real.(λ) .< 0.0)
    classification = is_stable ? "Asymptotically Stable (Sink)" : "Unstable (Saddle/Source)"
    
    println("\nFixed Point #$i:")
    println("  Position (x*, y*) : $(fp[1:2]) m")
    println("  Velocity (vx*, vy*): $(fp[3:4]) m/s")
    println("\n  Jacobian Matrix J (4x4):")
    display(round.(J, digits=4))
    
    println("\n  Eigenvalues (λ):")
    for (j, val) in enumerate(λ)
        re, im = round(real(val), digits=4), round(imag(val), digits=4)
        sign_str = im >= 0 ? "+" : "-"
        println("    λ_$(j) = $re $sign_str $(abs(im))i")
    end
        
    println("\n  Classification: $classification")
    println("-"^50)
end

# ---- Plotting Equilibrium Points & System Trajectory ---- #

fp_x = [fp[1] for fp in fixed_points]
fp_y = [fp[2] for fp in fixed_points]

p_summary = plot(
    x, y,
    label="Trajectory",
    lw=1.5,
    linestyle=:dash,
    alpha=0.7,
    color=:blue,
    xlabel="x (m)",
    ylabel="y (m)",
    title="Spring-Mass System: Trajectory & Clustered Fixed Points",
    aspect_ratio=:equal,
    legend=:topright
)

scatter!(
    p_summary,
    [params.rₐ[1], params.rᵦ[1]],
    [params.rₐ[2], params.rᵦ[2]],
    markersize=8,
    color=:black,
    markerstrokewidth=0,
    label="Anchors"
)

scatter!(
    p_summary,
    fp_x,
    fp_y,
    markersize=10,
    color=:red,
    shape=:star5,
    markerstrokewidth=1,
    label="Clustered Fixed Points"
)

for (i, fp) in enumerate(fixed_points)
    annotate!(p_summary, fp[1] + 0.3, fp[2] + 0.3, text("u*$i", :red, 9, :bold))
end

display(p_summary)
savefig("Problem_13/figures/fixed_points_analysis.png")
