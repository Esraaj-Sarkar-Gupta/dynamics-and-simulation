module Animations

using Plots

export create_animation, create_trajectory_plot, save_snapshot

"""
Generates an animated GIF of the triple pendulum and saves it to output_dir.
"""
function create_animation(sol, params, output_dir::String; filename="triple_pendulum.gif", fps=30)
    mkpath(output_dir)
    gif_path = joinpath(output_dir, filename)

    x1 = [u[1] for u in sol.u]; y1 = [u[2] for u in sol.u]
    x2 = [u[3] for u in sol.u]; y2 = [u[4] for u in sol.u]
    x3 = [u[5] for u in sol.u]; y3 = [u[6] for u in sol.u]

    lim = (params.l1 + params.l2 + params.l3) * 1.15

    anim = @animate for i in 1:length(sol.t)
        p = plot(
            xlims=(-lim, lim),
            ylims=(-lim, lim),
            aspect_ratio=:equal,
            legend=:topright,
            title="Triple Pendulum Simulation (t = $(round(sol.t[i], digits=2)) s)",
            xlabel="x (m)",
            ylabel="y (m)",
            grid=true
        )

        # Plot path traces up to current frame
        plot!(p, x1[1:i], y1[1:i], label="Mass 1 Path", color=:blue, alpha=0.35, lw=1.2)
        plot!(p, x2[1:i], y2[1:i], label="Mass 2 Path", color=:orange, alpha=0.35, lw=1.2)
        plot!(p, x3[1:i], y3[1:i], label="Mass 3 Path", color=:green, alpha=0.35, lw=1.2)

        # Pendulum rods
        px = [0.0, x1[i], x2[i], x3[i]]
        py = [0.0, y1[i], y2[i], y3[i]]
        plot!(p, px, py, color=:black, lw=2.0, label="")

        # Pivot point
        scatter!(p, [0.0], [0.0], color=:black, markersize=5, label="")

        # Masses
        scatter!(p, [x1[i]], [y1[i]], color=:blue, markersize=7, label="Mass 1")
        scatter!(p, [x2[i]], [y2[i]], color=:orange, markersize=7, label="Mass 2")
        scatter!(p, [x3[i]], [y3[i]], color=:green, markersize=7, label="Mass 3")
    end

    gif(anim, gif_path, fps=fps)
    println("Saved animation GIF to: ", gif_path)
    return gif_path
end

"""
Plots and saves the full trajectory paths for all 3 end effectors.
"""
function create_trajectory_plot(sol, params, output_dir::String; filename="trajectories.png")
    mkpath(output_dir)
    plot_path = joinpath(output_dir, filename)

    x1 = [u[1] for u in sol.u]; y1 = [u[2] for u in sol.u]
    x2 = [u[3] for u in sol.u]; y2 = [u[4] for u in sol.u]
    x3 = [u[5] for u in sol.u]; y3 = [u[6] for u in sol.u]

    lim = (params.l1 + params.l2 + params.l3) * 1.15

    p = plot(
        xlims=(-lim, lim),
        ylims=(-lim, lim),
        aspect_ratio=:equal,
        legend=:topright,
        title="Triple Pendulum End-Effector Trajectories",
        xlabel="x (m)",
        ylabel="y (m)",
        grid=true
    )

    plot!(p, x1, y1, label="Mass 1", color=:blue, lw=1.5)
    plot!(p, x2, y2, label="Mass 2", color=:orange, lw=1.5)
    plot!(p, x3, y3, label="Mass 3", color=:green, lw=1.5)
    scatter!(p, [0.0], [0.0], color=:black, markersize=6, label="Pivot")

    savefig(p, plot_path)
    println("Saved trajectory plot to: ", plot_path)
    return plot_path
end

"""
Saves a single frame snapshot from the simulation for inclusion in a report.
"""
function save_snapshot(sol, params, output_dir::String; frame_idx=nothing, filename="snapshot.png")
    mkpath(output_dir)
    snapshot_path = joinpath(output_dir, filename)

    if isnothing(frame_idx)
        frame_idx = max(1, div(length(sol.t), 2))
    end

    x1 = [u[1] for u in sol.u]; y1 = [u[2] for u in sol.u]
    x2 = [u[3] for u in sol.u]; y2 = [u[4] for u in sol.u]
    x3 = [u[5] for u in sol.u]; y3 = [u[6] for u in sol.u]

    lim = (params.l1 + params.l2 + params.l3) * 1.15
    i = frame_idx

    p = plot(
        xlims=(-lim, lim),
        ylims=(-lim, lim),
        aspect_ratio=:equal,
        legend=:topright,
        title="Triple Pendulum Snapshot (t = $(round(sol.t[i], digits=2)) s)",
        xlabel="x (m)",
        ylabel="y (m)",
        grid=true
    )

    plot!(p, x1[1:i], y1[1:i], label="Mass 1 Path", color=:blue, alpha=0.35, lw=1.2)
    plot!(p, x2[1:i], y2[1:i], label="Mass 2 Path", color=:orange, alpha=0.35, lw=1.2)
    plot!(p, x3[1:i], y3[1:i], label="Mass 3 Path", color=:green, alpha=0.35, lw=1.2)

    px = [0.0, x1[i], x2[i], x3[i]]
    py = [0.0, y1[i], y2[i], y3[i]]
    plot!(p, px, py, color=:black, lw=2.0, label="")

    scatter!(p, [0.0], [0.0], color=:black, markersize=5, label="Pivot")
    scatter!(p, [x1[i]], [y1[i]], color=:blue, markersize=7, label="Mass 1")
    scatter!(p, [x2[i]], [y2[i]], color=:orange, markersize=7, label="Mass 2")
    scatter!(p, [x3[i]], [y3[i]], color=:green, markersize=7, label="Mass 3")

    savefig(p, snapshot_path)
    println("Saved report snapshot to: ", snapshot_path)
    return snapshot_path
end

end # module Animations