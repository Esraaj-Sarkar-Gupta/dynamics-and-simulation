# ===== Double Pendulum: Visualization ===== #

module Animations

using Plots
using ..System: angles_to_cartesian, total_energy

export create_animation, create_trajectory_plot, save_snapshot, create_energy_plot

"""Extract the (x,y) trajectories of both point masses from an ODE solution."""
function extract_positions(sol, params)
    xy = [angles_to_cartesian(u[1], u[2], params) for u in sol.u]
    x1 = [p[1] for p in xy]; y1 = [p[2] for p in xy]
    x2 = [p[3] for p in xy]; y2 = [p[4] for p in xy]
    return x1, y1, x2, y2
end

"""Generate a GIF, including rod motion and the accumulated mass traces."""
function create_animation(sol, params, output_dir::String;
                          filename="double_pendulum.gif", fps=30)
    mkpath(output_dir)
    gif_path = joinpath(output_dir, filename)

    x1, y1, x2, y2 = extract_positions(sol, params)
    lim = 1.15 * (params.l1 + params.l2)

    anim = @animate for i in eachindex(sol.t)
        p = plot(
            xlims=(-lim, lim), ylims=(-lim, lim),
            aspect_ratio=:equal,
            legend=:topright,
            title="Double Pendulum -- Newton--Euler (t = $(round(sol.t[i], digits=2)) s)",
            xlabel="x (m)", ylabel="y (m)", grid=true
        )

        plot!(p, x1[1:i], y1[1:i], color=:blue, alpha=0.35,
              lw=1.2, label="Mass 1 Path")
        plot!(p, x2[1:i], y2[1:i], color=:orange, alpha=0.35,
              lw=1.2, label="Mass 2 Path")

        plot!(p, [0.0, x1[i], x2[i]], [0.0, y1[i], y2[i]],
              color=:black, lw=2, label="")
        scatter!(p, [0.0], [0.0], color=:black, markersize=5, label="")
        scatter!(p, [x1[i]], [y1[i]], color=:blue, markersize=7, label="Mass 1")
        scatter!(p, [x2[i]], [y2[i]], color=:orange, markersize=7, label="Mass 2")
    end

    gif(anim, gif_path, fps=fps)
    println("Saved animation GIF to: ", gif_path)
    return gif_path
end

"""Plot the complete Cartesian trajectories of both masses."""
function create_trajectory_plot(sol, params, output_dir::String;
                                filename="trajectories.png")
    mkpath(output_dir)
    plot_path = joinpath(output_dir, filename)

    x1, y1, x2, y2 = extract_positions(sol, params)
    lim = 1.15 * (params.l1 + params.l2)

    p = plot(
        xlims=(-lim, lim), ylims=(-lim, lim), aspect_ratio=:equal,
        legend=:topright, title="Double Pendulum End-Effector Trajectories",
        xlabel="x (m)", ylabel="y (m)", grid=true
    )
    plot!(p, x1, y1, color=:blue, lw=1.5, label="Mass 1")
    plot!(p, x2, y2, color=:orange, lw=1.5, label="Mass 2")
    scatter!(p, [0.0], [0.0], color=:black, markersize=6, label="Pivot")

    savefig(p, plot_path)
    println("Saved trajectory plot to: ", plot_path)
    return plot_path
end

"""Save a report-friendly snapshot with rods, masses and trajectory traces."""
function save_snapshot(sol, params, output_dir::String;
                       frame_idx=nothing, filename="report_snapshot.png")
    mkpath(output_dir)
    snapshot_path = joinpath(output_dir, filename)

    i = isnothing(frame_idx) ? max(1, div(length(sol.t), 2)) : frame_idx
    1 <= i <= length(sol.t) || throw(ArgumentError("frame_idx out of range"))

    x1, y1, x2, y2 = extract_positions(sol, params)
    lim = 1.15 * (params.l1 + params.l2)

    p = plot(
        xlims=(-lim, lim), ylims=(-lim, lim), aspect_ratio=:equal,
        legend=:topright,
        title="Double Pendulum Snapshot (t = $(round(sol.t[i], digits=2)) s)",
        xlabel="x (m)", ylabel="y (m)", grid=true
    )
    plot!(p, x1[1:i], y1[1:i], color=:blue, alpha=0.35,
          lw=1.2, label="Mass 1 Path")
    plot!(p, x2[1:i], y2[1:i], color=:orange, alpha=0.35,
          lw=1.2, label="Mass 2 Path")

    plot!(p, [0.0, x1[i], x2[i]], [0.0, y1[i], y2[i]],
          color=:black, lw=2.0, label="")
    scatter!(p, [0.0], [0.0], color=:black, markersize=5, label="Pivot")
    scatter!(p, [x1[i]], [y1[i]], color=:blue, markersize=7, label="Mass 1")
    scatter!(p, [x2[i]], [y2[i]], color=:orange, markersize=7, label="Mass 2")

    savefig(p, snapshot_path)
    println("Saved report snapshot to: ", snapshot_path)
    return snapshot_path
end

"""Plot relative energy drift to check the numerical integration."""
function create_energy_plot(sol, params, output_dir::String;
                            filename="energy_conservation.png")
    mkpath(output_dir)
    plot_path = joinpath(output_dir, filename)

    E = [total_energy(u, params) for u in sol.u]
    scale = max(abs(E[1]), 1.0)  # Avoid division by zero if E(0) = 0.
    drift = (E .- E[1]) ./ scale

    p = plot(sol.t, drift, lw=1.6, legend=false, grid=true,
             xlabel="Time (s)", ylabel="[E(t) - E(0)] / max(|E(0)|, 1 J)",
             title="Mechanical Energy Conservation")
    savefig(p, plot_path)
    println("Saved energy plot to: ", plot_path)
    return plot_path
end

end # module Animations
