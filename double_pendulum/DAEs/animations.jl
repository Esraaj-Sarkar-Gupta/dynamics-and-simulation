module Animations

using Plots

export plot_trajectory, animate_pendulum

function plot_trajectory(sol, save_dir::String)
    # Extract positions
    x1 = [u[1] for u in sol.u]
    y1 = [u[2] for u in sol.u]
    x2 = [u[3] for u in sol.u]
    y2 = [u[4] for u in sol.u]

    plt = plot(x1, y1, label="Mass 1", color=:blue, lw=1.5, aspect_ratio=:equal)
    plot!(plt, x2, y2, label="Mass 2", color=:red, lw=1.5)
    
    title!(plt, "Double Pendulum Trajectory (Cartesian DAE)")
    xlabel!(plt, "x (m)")
    ylabel!(plt, "y (m)")
    
    filepath = joinpath(save_dir, "trajectory.png")
    savefig(plt, filepath)
    println("Saved trajectory plot to $filepath")
end

function animate_pendulum(sol, p, save_dir::String; fps=30)
    t_span = sol.t[end] - sol.t[1]
    n_frames = round(Int, t_span * fps)
    t_frames = range(sol.t[1], sol.t[end], length=n_frames)
    
    l_max = p.l1 + p.l2
    limit = l_max * 1.1

    anim = @animate for t in t_frames
        state = sol(t)
        x1, y1, x2, y2 = state[1], state[2], state[3], state[4]
        
        plot([0, x1, x2], [0, y1, y2], 
            lw=3, 
            color=:black, 
            label=false, 
            aspect_ratio=:equal,
            xlims=(-limit, limit), 
            ylims=(-limit, limit)
        )
        
        scatter!([0], [0], color=:black, ms=6, label="Pivot")
        scatter!([x1], [y1], color=:blue, ms=8, label="m1")
        scatter!([x2], [y2], color=:red, ms=8, label="m2")
        title!("Double Pendulum t = $(round(t, digits=2))s")
    end

    filepath = joinpath(save_dir, "double_pendulum.mp4")
    mp4(anim, filepath, fps=fps)
    println("Saved animation to $filepath")
end

end # module