# ===== Triple Pendulum: Plots and Animations ===== #

module Animations

using Plots

export create_animation, create_trajectory_plot, create_angle_plot,
       save_snapshot, create_energy_plot

# Obtain Cartesian coordinates from the angular ODE state.
function cartesian_tracks(sol, p)
    x1 = Float64[]; y1 = Float64[]
    x2 = Float64[]; y2 = Float64[]
    x3 = Float64[]; y3 = Float64[]
    for u in sol.u
        θ1, θ2, θ3 = u[1:3]
        xx1 = p.l1*sin(θ1)
        yy1 = -p.l1*cos(θ1)
        xx2 = xx1 + p.l2*sin(θ2)
        yy2 = yy1 - p.l2*cos(θ2)
        xx3 = xx2 + p.l3*sin(θ3)
        yy3 = yy2 - p.l3*cos(θ3)
        push!(x1, xx1); push!(y1, yy1)
        push!(x2, xx2); push!(y2, yy2)
        push!(x3, xx3); push!(y3, yy3)
    end
    return (x1=x1, y1=y1, x2=x2, y2=y2, x3=x3, y3=y3)
end

function axes_base(p; title="Triple Pendulum -- Newton--Euler")
    lim = 1.15 * (p.l1 + p.l2 + p.l3)
    return plot(xlims=(-lim, lim), ylims=(-lim, lim), aspect_ratio=:equal,
                xlabel="x (m)", ylabel="y (m)", title=title,
                legend=:topright, grid=true)
end

function draw_frame!(fig, tracks, i; trails=true)
    if trails
        plot!(fig, tracks.x1[1:i], tracks.y1[1:i],
              color=:blue, alpha=0.40, lw=1.2, label="Mass 1 path")
        plot!(fig, tracks.x2[1:i], tracks.y2[1:i],
              color=:orange, alpha=0.40, lw=1.2, label="Mass 2 path")
        plot!(fig, tracks.x3[1:i], tracks.y3[1:i],
              color=:green, alpha=0.40, lw=1.2, label="Mass 3 path")
    end
    plot!(fig, [0.0, tracks.x1[i], tracks.x2[i], tracks.x3[i]],
               [0.0, tracks.y1[i], tracks.y2[i], tracks.y3[i]],
          color=:black, lw=2.5, label="")
    scatter!(fig, [0.0], [0.0], color=:black, markersize=5, label="")
    scatter!(fig, [tracks.x1[i]], [tracks.y1[i]], color=:blue,
             markersize=7, label="Mass 1")
    scatter!(fig, [tracks.x2[i]], [tracks.y2[i]], color=:orange,
             markersize=7, label="Mass 2")
    scatter!(fig, [tracks.x3[i]], [tracks.y3[i]], color=:green,
             markersize=7, label="Mass 3")
    return fig
end

"""Plot and save the complete spatial trajectories of all three masses."""
function create_trajectory_plot(sol, p, output_dir::String;
                                filename="trajectories.png")
    mkpath(output_dir)
    tr = cartesian_tracks(sol, p)
    fig = axes_base(p; title="Triple Pendulum End-Effector Trajectories (Newton--Euler)")
    plot!(fig, tr.x1, tr.y1, color=:blue, lw=1.5, label="Mass 1")
    plot!(fig, tr.x2, tr.y2, color=:orange, lw=1.5, label="Mass 2")
    plot!(fig, tr.x3, tr.y3, color=:green, lw=1.5, label="Mass 3")
    scatter!(fig, [0.0], [0.0], color=:black, markersize=5, label="Pivot")
    path = joinpath(output_dir, filename)
    savefig(fig, path)
    println("Saved trajectories: ", path)
    return path
end

"""Plot θ1(t), θ2(t), θ3(t) in radians."""
function create_angle_plot(sol, p, output_dir::String;
                           filename="angles_plot.png")
    mkpath(output_dir)
    fig = plot(sol.t, [u[1] for u in sol.u], label="θ₁", lw=1.5,
               xlabel="Time (s)", ylabel="Angle (rad)", grid=true,
               title="Triple Pendulum Angular Motion (Newton--Euler)")
    plot!(fig, sol.t, [u[2] for u in sol.u], label="θ₂", lw=1.5)
    plot!(fig, sol.t, [u[3] for u in sol.u], label="θ₃", lw=1.5)
    path = joinpath(output_dir, filename)
    savefig(fig, path)
    println("Saved angles: ", path)
    return path
end

"""Save a report-friendly frame, including trajectories up to that time."""
function save_snapshot(sol, p, output_dir::String;
                       frame_idx=nothing, filename="report_snapshot.png")
    mkpath(output_dir)
    i = isnothing(frame_idx) ? max(1, div(length(sol.t), 2)) :
        clamp(frame_idx, 1, length(sol.t))
    tr = cartesian_tracks(sol, p)
    fig = axes_base(p; title="Triple Pendulum Snapshot (t=$(round(sol.t[i], digits=2)) s)")
    draw_frame!(fig, tr, i)
    path = joinpath(output_dir, filename)
    savefig(fig, path)
    println("Saved snapshot: ", path)
    return path
end

"""Create an animated GIF; frame_step=2 keeps the output reasonably small."""
function create_animation(sol, p, output_dir::String;
                          filename="triple_pendulum.gif", fps=30, frame_step=2)
    mkpath(output_dir)
    frame_step >= 1 || throw(ArgumentError("frame_step must be positive"))
    tr = cartesian_tracks(sol, p)
    anim = @animate for i in 1:frame_step:length(sol.t)
        fig = axes_base(p; title="Triple Pendulum (t=$(round(sol.t[i], digits=2)) s)")
        draw_frame!(fig, tr, i)
    end
    path = joinpath(output_dir, filename)
    gif(anim, path, fps=fps)
    println("Saved GIF: ", path)
    return path
end

"""Plot relative energy drift; should stay near zero for a converged solution."""
function create_energy_plot(sol, p, output_dir::String;
                            filename="energy_conservation.png")
    mkpath(output_dir)
    # System module is loaded by main.jl; this module avoids depending on it.
    # Energy computed directly from the Cartesian velocities of the point masses.
    function energy(u)
        θ = u[1:3]; ω = u[4:6]
        l = (p.l1, p.l2, p.l3)
        m = (p.m1, p.m2, p.m3)
        vx, vy, y, e = 0.0, 0.0, 0.0, 0.0
        for j in 1:3
            vx += l[j]*cos(θ[j])*ω[j]
            vy += l[j]*sin(θ[j])*ω[j]
            y  -= l[j]*cos(θ[j])
            e += 0.5*m[j]*(vx^2+vy^2) + m[j]*p.g*y
        end
        return e
    end
    E = [energy(u) for u in sol.u]
    scale = max(abs(E[1]), 1.0) # avoids division by zero at E(0)=0
    drift = (E .- E[1]) ./ scale
    fig = plot(sol.t, drift, lw=1.5, legend=false, xlabel="Time (s)",
               ylabel="[E(t)-E(0)] / max(|E(0)|,1 J)", grid=true,
               title="Mechanical Energy Drift (Newton--Euler)")
    path = joinpath(output_dir, filename)
    mkpath(output_dir)
    savefig(fig, path)
    println("Saved energy check: ", path)
    println("Maximum absolute scaled energy drift: ", maximum(abs.(drift)))
    return path
end

end # module Animations
