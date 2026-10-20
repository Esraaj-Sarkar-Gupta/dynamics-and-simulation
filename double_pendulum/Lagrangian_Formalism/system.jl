module System

export DoublePenduParams, double_pendulum

Base.@kwdef struct DoublePenduParams
    m1 :: Float64 = 1.0
    m2 :: Float64 = 1.0
    l1 :: Float64 = 1.0
    l2 :: Float64 = 1.0
    g  :: Float64 = 9.810
end

function double_pendulum(du, u, p::DoublePenduParams, t)
    θ1, θ2, ω1, ω2 = u
    m1, m2, l1, l2, g = p.m1, p.m2, p.l1, p.l2, p.g

    dθ1 = ω1
    dθ2 = ω2

    dω1 = (g*m1*sin(θ1) + g*m2*sin(θ1) - cos(θ1 - θ2)*g*m2*sin(θ2) +
            (dθ2^2)*l2*m2*sin(θ1 - θ2) +
            (dθ1^2)*cos(θ1 - θ2)*l1*m2*sin(θ1 - θ2)) /
            (-l1*m1 - l1*m2 + (cos(θ1 - θ2)^2)*l1*m2)

    dω2 = (g*m1*sin(θ2) + g*m2*sin(θ2) - cos(θ1 - θ2)*g*m1*sin(θ1) - cos(θ1 - θ2)*g*m2*sin(θ1) -
            (dθ1^2)*l1*m1*sin(θ1 - θ2) - (dθ1^2)*l1*m2*sin(θ1 - θ2) -
            (dθ2^2)*cos(θ1 - θ2)*l2*m2*sin(θ1 - θ2)) / 
            (-l2*m1 - l2*m2 + (cos(θ1 - θ2)^2)*l2*m2)

    # In-place updates
    du[1] = dθ1
    du[2] = dθ2
    du[3] = dω1
    du[4] = dω2

    return nothing
end

end # Module