module FourBarDAESystem

using LinearAlgebra
using DifferentialEquations
using ..FourBarDAESymbolics: derive_cartesian_dae

export FourBarParameters, FourBarModel, build_system, initial_state,
       constraint_residual, project_state!, kkt_system, dynamics, simulate,
       configuration, total_energy

Base.@kwdef struct FourBarParameters
    l0::Float64 = 1.5
    l1::Float64 = 0.6
    l2::Float64 = 1.4
    l3::Float64 = 1.0
    m1::Float64 = 1.0
    m2::Float64 = 1.2
    m3::Float64 = 0.9
    g::Float64 = 9.81
end

parameter_vector(p::FourBarParameters) =
    [p.l0, p.l1, p.l2, p.l3, p.m1, p.m2, p.m3, p.g]

struct FourBarModel{S,F}
    parameters::FourBarParameters
    symbolic::S
    parameter_values::Vector{Float64}
    torque::F
end

function build_system(p::FourBarParameters=FourBarParameters(); torque=t -> 0.0)
    all(x -> isfinite(x) && x > 0,
        (p.l0, p.l1, p.l2, p.l3, p.m1, p.m2, p.m3)) ||
        throw(ArgumentError("lengths and masses must be positive and finite"))
    isfinite(p.g) && p.g >= 0 ||
        throw(ArgumentError("gravity must be finite and nonnegative"))
    return FourBarModel(p, derive_cartesian_dae(), parameter_vector(p), torque)
end

"""Position constraints Ψ(q), in m², in the report's link order."""
function constraint_residual(model::FourBarModel, q::AbstractVector)
    return Vector{Float64}(model.symbolic.functions.Psi(
        q, zeros(4), model.parameter_values))
end

"""Find a closed configuration from input angle and circle-intersection branch."""
function initial_state(model::FourBarModel; theta1::Real=pi/3,
                       omega1::Real=0.0, branch::Int=1)
    branch in (-1, 1) || throw(ArgumentError("branch must be +1 or -1"))
    p = model.parameters
    rA = p.l1 .* [cos(theta1), sin(theta1)]
    dx, dy = p.l0-rA[1], -rA[2]
    alpha = atan(dy, dx)

    phi(delta) = begin
        theta2 = alpha + branch*delta
        xB = rA[1] + p.l2*cos(theta2)
        yB = rA[2] + p.l2*sin(theta2)
        ((xB-p.l0)^2 + yB^2 - p.l3^2)/2
    end
    flo, fhi = phi(0.0), phi(pi)
    scale = max(p.l0^2, p.l1^2, p.l2^2, p.l3^2)
    flo <= 1e-12*scale && fhi >= -1e-12*scale ||
        throw(DomainError(theta1, "four-bar cannot close at this input angle"))
    lo, hi = 0.0, pi
    for _ in 1:55
        mid = (lo+hi)/2
        if phi(mid) > 0
            hi = mid
        else
            lo = mid
        end
    end
    theta2 = alpha + branch*(lo+hi)/2
    rB = rA + p.l2 .* [cos(theta2), sin(theta2)]

    vA = p.l1*omega1 .* [-sin(theta1), cos(theta1)]
    dAB, dDB = rB-rA, rB-[p.l0, 0.0]
    # The two remaining velocity constraints determine vB at a regular pose.
    C = [dAB[1] dAB[2]; dDB[1] dDB[2]]
    rhs = [dot(dAB, vA), 0.0]
    abs(det(C)) > 1e-12*max(p.l2*p.l3, 1.0) ||
        throw(DomainError(theta1, "initial four-bar pose is at a toggle"))
    vB = C \ rhs
    return Float64[rA[1], rA[2], rB[1], rB[2], vA[1], vA[2], vB[1], vB[2]]
end

"""Project positions to Ψ=0 and velocities to J(q)q̇=0."""
function project_state!(u::AbstractVector, model::FourBarModel;
                        maxiter::Int=12, tol::Float64=1e-12)
    f, pars = model.symbolic.functions, model.parameter_values
    q, qd = @view(u[1:4]), @view(u[5:8])
    length_scale = max(model.parameters.l0, model.parameters.l1,
                       model.parameters.l2, model.parameters.l3)
    threshold = tol*length_scale^2
    for _ in 1:maxiter
        psi = Vector{Float64}(f.Psi(q, qd, pars))
        norm(psi, Inf) <= threshold && break
        J = Matrix{Float64}(f.J(q, qd, pars))
        gram = J*J'
        minimum(svdvals(gram)) > 1e-24*length_scale^4 ||
            throw(DomainError(q, "constraint Jacobian is rank deficient"))
        q .-= J' * (gram \ psi)
    end
    norm(constraint_residual(model, q), Inf) <= threshold ||
        throw(ErrorException("position constraint projection failed"))
    J = Matrix{Float64}(f.J(q, qd, pars))
    gram = J*J'
    qd .-= J' * (gram \ (J*qd))
    return u
end

"""Return positions, link angles, and angular rates for plotting/output."""
function configuration(model::FourBarModel, u::AbstractVector)
    p = model.parameters
    q, qd = u[1:4], u[5:8]
    rA, rB = q[1:2], q[3:4]
    vA, vB = qd[1:2], qd[3:4]
    dAB, dDB = rB-rA, rB-[p.l0, 0.0]
    cross2(a,b) = a[1]*b[2]-a[2]*b[1]
    angles = [atan(rA[2], rA[1]), atan(dAB[2], dAB[1]), atan(dDB[2], dDB[1])]
    rates = [cross2(rA,vA)/p.l1^2, cross2(dAB,vB-vA)/p.l2^2,
             cross2(dDB,vB)/p.l3^2]
    return (; rA, rB, rD=[p.l0,0.0], angles, rates)
end

"""Assemble K=[M_c J_c'; J_c 0] and its right-hand side."""
function kkt_system(model::FourBarModel, u::AbstractVector, t::Real)
    q, qd = u[1:4], u[5:8]
    f, pars = model.symbolic.functions, model.parameter_values # Symbolic Lib Helpers

    M = Matrix{Float64}(f.M(q, qd, pars))
    c = Vector{Float64}(f.c(q, qd, pars))
    J = Matrix{Float64}(f.J(q, qd, pars))
    Jdot = Matrix{Float64}(f.Jdot(q, qd, pars))
    gamma = Jdot*qd # J̇(q,q̇)q̇
    tau = Float64(model.torque(t))
    Q = (tau/model.parameters.l1^2) .* [-q[2], q[1], 0.0, 0.0]

    # Our big matrix
    K = zeros(7, 7)
    K[1:4, 1:4] .= M
    K[1:4, 5:7] .= J'
    K[5:7, 1:4] .= J
    rhs = [Q-c; -gamma]
    return K, rhs
end

function dynamics(model::FourBarModel, u::AbstractVector, t::Real)
    K, rhs = kkt_system(model, u, t)
    solution = K \ rhs
    return (; acceleration=solution[1:4], lambda=solution[5:7], matrix=K, rhs)
end

function rhs!(du, u, model::FourBarModel, t)
    du[1:4] .= u[5:8]
    du[5:8] .= dynamics(model, u, t).acceleration
    return nothing
end

"""Integrate the index-reduced constrained equations and project accepted states."""
function simulate(model::FourBarModel;
                  theta1::Real=pi/3, omega1::Real=0.0, branch::Int=1,
                  u0::Union{Nothing,AbstractVector}=nothing,
                  tspan::Tuple{<:Real,<:Real}=(0.0, 8.0), fps::Int=30,
                  reltol::Real=1e-9, abstol::Real=1e-11)
    fps > 0 || throw(ArgumentError("fps must be positive"))
    tspan[2] > tspan[1] || throw(ArgumentError("tspan must increase"))
    state = isnothing(u0) ? initial_state(model; theta1, omega1, branch) : Float64.(u0)
    length(state) == 8 || throw(DimensionMismatch("u0 must have 8 entries: [q; q̇]"))
    project_state!(state, model)
    prob = ODEProblem(rhs!, state, (Float64(tspan[1]), Float64(tspan[2])), model)
    callback = DiscreteCallback((u,t,integrator)->true,
        integrator->project_state!(integrator.u, model); save_positions=(false,false))
    n = max(2, round(Int, (tspan[2]-tspan[1])*fps)+1)
    sampled = range(Float64(tspan[1]), Float64(tspan[2]); length=n)
    raw = solve(prob, Tsit5(); callback, saveat=sampled, save_everystep=false,
                reltol, abstol)
    states = [project_state!(copy(x), model) for x in raw.u]
    return (; t=Float64.(raw.t), u=states, retcode=raw.retcode)
end

total_energy(model::FourBarModel, u::AbstractVector) =
    Float64(model.symbolic.functions.energy(u[1:4], u[5:8], model.parameter_values))

end # module FourBarDAESystem
