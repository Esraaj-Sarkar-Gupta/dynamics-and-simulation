module Animate

using Plots
using DifferentialEquations
using ..System

export plot_trajectories, plot_phase_space, plot_cartesian_path, render_snapshot, animate_pendulum

# Forward Kinematics
x1(θ1, l1) = l1 * sin(θ1)
y1(θ1, l1) = -l1 * cos(θ1)
x2(θ1, θ2, l1, l2) = x1(θ1, l1) + l2 * sin(θ2)
y2(θ1, θ2, l1, l2) = y1(θ1, l1) - l2 * cos(θ2)

"""
    plot_trajectories(sol, output_dir)

Plots angular position trajectories over time for θ₁ and θ₂.
"""
function plot_trajectories(sol, output_dir::String)
    plt = plot(
        sol.t, [sol[1, :], sol[2, :]],
        label=["θ₁ (rad)" "θ₂ (rad)"],
        xlabel="Time t (s)",
        ylabel="Angle (rad)",
        title="Double Pendulum Angular Trajectories",
        lw=1.5,
        legend=:outerright
    )
    savefig(plt, joinpath(output_dir, "angular_trajectories.png"))
    return plt
end

"""
    plot_phase_space(sol, output_dir)

Plots phase space diagrams (angle vs. angular velocity) for both links.
"""
function plot_phase_space(sol, output_dir::String)
    plt = plot(
        sol[1, :], sol[3, :],
        label="Link 1 (θ₁ vs ω₁)",
        xlabel="Angle θ (rad)",
        ylabel="Angular Velocity ω (rad/s)",
        title="Phase Space Dynamics",
        lw=1.2
    )
    plot!(plt, sol[2, :], sol[4, :], label="Link 2 (θ₂ vs ω₂)", lw=1.2)
    savefig(plt, joinpath(output_dir, "phase_space.png"))
    return plt
end

"""
    plot_cartesian_path(sol, params, output_dir)

Plots the 2D Cartesian spatial path executed by the end-effector mass m₂.
"""
function plot_cartesian_path(sol, params::DoublePenduParams, output_dir::String)
    θ1_vals = sol[1, :]
    θ2_vals = sol[2, :]

    px2 = x2.(θ1_vals, θ2_vals, params.l1, params.l2)
    py2 = y2.(θ1_vals, θ2_vals, params.l1, params.l2)

    plt = plot(
        px2, py2,
        label="Tip Trajectory",
        xlabel="x (m)",
        ylabel="y (m)",
        title="End-Effector Cartesian Path",
        aspect_ratio=:equal,
        lw=1.0,
        color=:crimson
    )
    scatter!(plt, [0], [0], color=:black, ms=5, label="Pivot")
    savefig(plt, joinpath(output_dir, "end_effector_trajectory.png"))
    return plt
end

"""
    render_frame(t, sol, params; trail_length=100, t_steps=nothing)

Internal helper function to render a single physical frame at time `t`.
"""
function render_frame(t::Float64, sol, params::DoublePenduParams; trail_length=100, t_steps=nothing, i=nothing)
    u = sol(t)
    θ1, θ2 = u[1], u[2]

    # Joint Positions
    px1, py1 = x1(θ1, params.l1), y1(θ1, params.l1)
    px2, py2 = x2(θ1, θ2, params.l1, params.l2), y2(θ1, θ2, params.l1, params.l2)

    # End-effector trail calculation
    if !isnothing(t_steps) && !isnothing(i)
        idx_start = max(1, i - trail_length)
        trail_t = t_steps[idx_start:i]
    else
        t_start = max(sol.t[1], t - 2.0)
        trail_t = range(t_start, t, length=50)
    end
    
    trail_u = [sol(tk) for tk in trail_t]
    tx2 = [x2(uk[1], uk[2], params.l1, params.l2) for uk in trail_u]
    ty2 = [y2(uk[1], uk[2], params.l1, params.l2) for uk in trail_u]

    lim = params.l1 + params.l2 + 0.2

    plt = plot(
        tx2, ty2,
        color=:lightskyblue,
        lw=1.5,
        label=nothing,
        xlims=(-lim, lim),
        ylims=(-lim, lim),
        aspect_ratio=:equal,
        title="Double Pendulum Dynamics (t = $(round(t, digits=2))s)",
        xlabel="x (m)",
        ylabel="y (m)",
        grid=true
    )

    # Link structure
    plot!(plt, [0, px1, px2], [0, py1, py2], lw=3, color=:black, label=nothing)

    # Mass joints
    scatter!(plt, [0], [0], color=:black, ms=6, label=nothing)
    scatter!(plt, [px1], [py1], color=:firebrick, ms=8, label="m₁")
    scatter!(plt, [px2], [py2], color=:darkorange, ms=8, label="m₂")

    return plt
end

"""
    render_snapshot(sol, params, output_dir; t_snap=sol.t[end]*0.35)

Saves a static image frame of the physical double pendulum state.
"""
function render_snapshot(sol, params::DoublePenduParams, output_dir::String; t_snap::Float64=sol.t[end]*0.35)
    plt = render_frame(t_snap, sol, params)
    savefig(plt, joinpath(output_dir, "animation_screenshot.png"))
    return plt
end

"""
    animate_pendulum(sol, params, output_dir; fps=30, trail_length=100)

Generates and saves a full animated GIF of the motion.
"""
function animate_pendulum(sol, params::DoublePenduParams, output_dir::String; fps::Int=30, trail_length::Int=100)
    tspan = (sol.t[1], sol.t[end])
    t_steps = range(tspan[1], tspan[2], step=1/fps)

    anim = @animate for i in 1:length(t_steps)
        t = t_steps[i]
        render_frame(t, sol, params; trail_length=trail_length, t_steps=t_steps, i=i)
    end

    gif_path = joinpath(output_dir, "double_pendulum_animation.gif")
    gif(anim, gif_path, fps=fps)
    return gif_path
end

end # Module