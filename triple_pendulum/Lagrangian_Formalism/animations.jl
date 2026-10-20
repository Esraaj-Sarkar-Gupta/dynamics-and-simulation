# animations.jl
# Time-series plots, spatial trajectory tracking, and static report snapshots

using Plots

function plot_and_animate(sol, p, output_dir::String; snapshot_time::Float64=2.5)
    println("5. Generating plots, spatial trajectories, and static report graphics...")
    
    # Ensure child directory exists
    mkpath(output_dir)

    l1, l2, l3 = p[4], p[5], p[6]
    l_total = l1 + l2 + l3

    # ---- 1. Angle Trajectory Plot ---- #
    p_traj = plot(sol, vars=[1, 2, 3], 
                  label=["θ₁" "θ₂" "θ₃"], 
                  linewidth=1.5, 
                  title="Triple Pendulum Angular Displacement",
                  xlabel="Time (s)", 
                  ylabel="Angle (rad)",
                  legend=:bottomleft,
                  gridalpha=0.4)

    plot_path = joinpath(output_dir, "angles_plot.png")
    savefig(p_traj, plot_path)
    println("   ✓ Saved angle time-series plot to: ", plot_path)

    # ---- 2. Compute Forward Kinematic Trajectories ---- #
    t_eval = range(sol.t[1], sol.t[end], length=1000)
    
    x1_path = Float64[]
    y1_path = Float64[]
    x2_path = Float64[]
    y2_path = Float64[]
    x3_path = Float64[]
    y3_path = Float64[]

    for t in t_eval
        u = sol(t)
        θ1, θ2, θ3 = u[1], u[2], u[3]
        
        # Forward Kinematics (zero angle hanging downwards)
        x1 = l1 * sin(θ1)
        y1 = -l1 * cos(θ1)
        
        x2 = x1 + l2 * sin(θ2)
        y2 = y1 - l2 * cos(θ2)
        
        x3 = x2 + l3 * sin(θ3)
        y3 = y2 - l3 * cos(θ3)

        push!(x1_path, x1); push!(y1_path, y1)
        push!(x2_path, x2); push!(y2_path, y2)
        push!(x3_path, x3); push!(y3_path, y3)
    end

    # ---- 3. Effector Trajectories Plot (xy-Plane) ---- #
    p_eff = plot(x1_path, y1_path, 
                 label="Mass 1 Path (m₁)", 
                 linewidth=1.2, 
                 alpha=0.7, 
                 color=:crimson)
    
    plot!(p_eff, x2_path, y2_path, 
          label="Mass 2 Path (m₂)", 
          linewidth=1.2, 
          alpha=0.7, 
          color=:darkorange)
    
    plot!(p_eff, x3_path, y3_path, 
          label="Mass 3 Path (m₃ - End Effector)", 
          linewidth=1.5, 
          alpha=0.9, 
          color=:royalblue)

    # Add origin anchor
    scatter!(p_eff, [0.0], [0.0], marker=:square, markersize=6, color=:black, label="Pivot Point")

    plot!(p_eff,
          title="Spatial Trajectories of All Three Effectors",
          xlabel="x Position (m)",
          ylabel="y Position (m)",
          xlims=(-l_total*1.15, l_total*1.15),
          ylims=(-l_total*1.15, l_total*1.15),
          aspect_ratio=:equal,
          legend=:topright,
          gridalpha=0.4)

    eff_path = joinpath(output_dir, "effector_trajectories.png")
    savefig(p_eff, eff_path)
    println("   Saved end-effector spatial trajectories to: ", eff_path)

    # ---- 4. Single Snapshot for Report ---- #
    u_snap = sol(snapshot_time)
    θ1, θ2, θ3 = u_snap[1], u_snap[2], u_snap[3]

    x1_s = l1 * sin(θ1)
    y1_s = -l1 * cos(θ1)
    x2_s = x1_s + l2 * sin(θ2)
    y2_s = y1_s - l2 * cos(θ2)
    x3_s = x2_s + l3 * sin(θ3)
    y3_s = y2_s - l3 * cos(θ3)

    x_nodes = [0.0, x1_s, x2_s, x3_s]
    y_nodes = [0.0, y1_s, y2_s, y3_s]

    p_snap = plot(x3_path[t_eval .<= snapshot_time], y3_path[t_eval .<= snapshot_time],
                  linewidth=1.0, 
                  linestyle=:dash, 
                  color=:gray)

    plot!(p_snap, x_nodes, y_nodes, 
          marker=:circle, 
          markersize=8, 
          linewidth=3, 
          color=:navy, 
          xlims=(-l_total*1.15, l_total*1.15),
          ylims=(-l_total*1.15, l_total*1.15),
          aspect_ratio=:equal,
          title="Physical State Snapshot (t = $(snapshot_time) s)",
          xlabel="x (m)",
          ylabel="y (m)",
          legend=:none,
          gridalpha=0.3)

    snap_path = joinpath(output_dir, "pendulum_snapshot.png")
    savefig(p_snap, snap_path)
    println("   Saved static report snapshot to: ", snap_path)

    # ---- 5. Physical Animation GIF ---- #
    fps = 30
    t_steps = range(sol.t[1], sol.t[end], length=Int(ceil((sol.t[end] - sol.t[1]) * fps)))
    
    anim = @animate for t in t_steps
        u = sol(t)
        θ1_i, θ2_i, θ3_i = u[1], u[2], u[3]
        
        x1_i = l1 * sin(θ1_i)
        y1_i = -l1 * cos(θ1_i)
        x2_i = x1_i + l2 * sin(θ2_i)
        y2_i = y1_i - l2 * cos(θ2_i)
        x3_i = x2_i + l3 * sin(θ3_i)
        y3_i = y2_i - l3 * cos(θ3_i)
        
        plot([0.0, x1_i, x2_i, x3_i], [0.0, y1_i, y2_i, y3_i], 
             marker=:circle, 
             markersize=6, 
             linewidth=3, 
             color=:navy,
             legend=false,
             xlims=(-l_total*1.15, l_total*1.15),
             ylims=(-l_total*1.15, l_total*1.15),
             aspect_ratio=:equal,
             title="Triple Pendulum Motion (t = $(round(t, digits=2)) s)",
             xlabel="x (m)",
             ylabel="y (m)",
             gridalpha=0.3)
    end

    anim_path = joinpath(output_dir, "triple_pendulum.gif")
    gif(anim, anim_path, fps=fps)
    println("   Saved physical animation to: ", anim_path)
end