module FourBarDAEAnimations

using CairoMakie
using Printf
using ..FourBarDAESystem: FourBarModel, configuration, total_energy

export save_visualizations

function unwrap_angles(x::AbstractVector)
    y = collect(Float64, x)
    for i in 2:length(y)
        while y[i]-y[i-1] > pi; y[i] -= 2pi; end
        while y[i]-y[i-1] < -pi; y[i] += 2pi; end
    end
    return y
end

"""Save joint trajectories, a midpoint snapshot, and an MP4 animation."""
function save_visualizations(model::FourBarModel, sol;
                             outdir::AbstractString="output", fps::Int=30)
    fps > 0 || throw(ArgumentError("fps must be positive"))
    mkpath(outdir)
    n = length(sol.t)
    n >= 2 || throw(ArgumentError("need at least two frames"))
    t = Float64.(sol.t)
    xA, yA, xB, yB = [zeros(n) for _ in 1:4]
    angles, rates, energy = zeros(3,n), zeros(3,n), zeros(n)
    for (i,u) in enumerate(sol.u)
        s = configuration(model, u)
        xA[i], yA[i] = s.rA
        xB[i], yB[i] = s.rB
        angles[:,i] .= s.angles
        rates[:,i] .= s.rates
        energy[i] = total_energy(model,u)
    end

    fig = Figure(size=(1200,840), fontsize=16)
    paths = Axis(fig[1:2,1]; title="Joint trajectories", xlabel="x (m)",
                 ylabel="y (m)", aspect=DataAspect())
    ang = Axis(fig[1,2]; title="Link orientations", xlabel="Time (s)", ylabel="Angle (rad)")
    rate = Axis(fig[2,2]; title="Angular velocities", xlabel="Time (s)", ylabel="rad/s")
    en = Axis(fig[3,1:2]; title="Mechanical energy", xlabel="Time (s)", ylabel="Energy (J)")
    lines!(paths,xA,yA;color=:royalblue,linewidth=2.5,label="A")
    lines!(paths,xB,yB;color=:darkorange,linewidth=2.5,label="B")
    scatter!(paths,[0.0,model.parameters.l0],[0.0,0.0];color=:black,markersize=11)
    axislegend(paths;position=:rb)
    colors = (:royalblue,:darkorange,:seagreen)
    for i in 1:3
        lines!(ang,t,unwrap_angles(view(angles,i,:));color=colors[i],label="θ$i",linewidth=2)
        lines!(rate,t,view(rates,i,:);color=colors[i],label="ω$i",linewidth=2)
    end
    axislegend(ang;position=:rb); axislegend(rate;position=:rb)
    lines!(en,t,energy;color=:purple,linewidth=2.3)
    trajectories_path = joinpath(outdir,"trajectories.png")
    save(trajectories_path,fig)

    O, D = Point2f(0,0), Point2f(model.parameters.l0,0)
    getA(i)=Point2f(xA[i],yA[i]); getB(i)=Point2f(xB[i],yB[i])
    link1=Observable(Point2f[O,getA(1)])
    link2=Observable(Point2f[getA(1),getB(1)])
    link3=Observable(Point2f[D,getB(1)])
    joints=Observable(Point2f[getA(1),getB(1)])
    traceA=Observable(Point2f[getA(1)]); traceB=Observable(Point2f[getB(1)])
    animfig=Figure(size=(900,700),fontsize=19)
    ax=Axis(animfig[1,1];aspect=DataAspect(),xlabel="x (m)",ylabel="y (m)",
            title="Four-bar linkage | Cartesian DAE")
    lines!(ax,traceA;color=(:royalblue,0.38),linewidth=2)
    lines!(ax,traceB;color=(:darkorange,0.38),linewidth=2)
    lines!(ax,[O,D];color=:black,linewidth=7)
    lines!(ax,link1;color=:royalblue,linewidth=9)
    lines!(ax,link2;color=:darkorange,linewidth=9)
    lines!(ax,link3;color=:seagreen,linewidth=9)
    scatter!(ax,[O,D];color=:black,markersize=19)
    scatter!(ax,joints;color=:black,markersize=16)
    allx=vcat([0.0,model.parameters.l0],xA,xB); ally=vcat([0.0,0.0],yA,yB)
    span=max(maximum(allx)-minimum(allx),maximum(ally)-minimum(ally),0.1)
    pad=max(0.15,0.12span)
    limits!(ax,minimum(allx)-pad,maximum(allx)+pad,minimum(ally)-pad,maximum(ally)+pad)
    function update!(i)
        A,B=getA(i),getB(i)
        link1[]=Point2f[O,A]; link2[]=Point2f[A,B]; link3[]=Point2f[D,B]
        joints[]=Point2f[A,B]
        traceA[]=Point2f[getA(k) for k in 1:i]
        traceB[]=Point2f[getB(k) for k in 1:i]
        ax.title=@sprintf("Four-bar linkage | t = %.2f s",t[i])
    end
    snapshot_path=joinpath(outdir,"snapshot.png")
    update!(cld(n,2)); save(snapshot_path,animfig)
    animation_path=joinpath(outdir,"fourbar_dae.mp4")
    record(animfig,animation_path,1:n;framerate=fps) do i; update!(i); end
    return (;trajectories_path,snapshot_path,animation_path,
            energy_min=minimum(energy),energy_max=maximum(energy))
end

end # module FourBarDAEAnimations
