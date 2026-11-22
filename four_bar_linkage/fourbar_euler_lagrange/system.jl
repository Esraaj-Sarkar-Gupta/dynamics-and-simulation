module FourBarSystem

using LinearAlgebra
using DifferentialEquations
using ..FourBarSymbolics: derive_euler_lagrange

export FourBarParameters, FourBarModel, build_system, body_positions,
       initial_angles, configuration, closure_residual, project_state!,
       dynamics, simulate, total_energy, constraint_multipliers

Base.@kwdef struct FourBarParameters
    l0::Float64 = 1.5   # O -> D, stationary
    l1::Float64 = 0.6   # O -> A
    l2::Float64 = 1.4   # A -> B
    l3::Float64 = 1.0   # D -> B
    m1::Float64 = 1.0
    m2::Float64 = 1.2
    m3::Float64 = 0.9
    g::Float64  = 9.81
end

parameter_vector(p::FourBarParameters) =
    [p.l0, p.l1, p.l2, p.l3, p.m1, p.m2, p.m3, p.g]

struct FourBarModel{S,F}
    parameters::FourBarParameters
    symbolic::S
    parameter_values::Vector{Float64}
    torque::F                  # torque on OA at O, tau(t)
end

"""Derive/compile symbolic expressions once for a set of parameters."""
function build_system(p::FourBarParameters=FourBarParameters();
                      torque=t -> 0.0)

    all(x -> isfinite(x) && x > 0,
        (p.l0, p.l1, p.l2, p.l3, p.m1, p.m2, p.m3)) ||
        throw(ArgumentError("lengths and masses must be positive and finite"))
    (isfinite(p.g) && p.g >= 0) ||
        throw(ArgumentError("gravity must be finite and nonnegative"))

    symbolic = derive_euler_lagrange()
    return FourBarModel(p, symbolic, parameter_vector(p), torque)
end

"""Forward kinematics; B is obtained directly from OA and AB."""
function body_positions(p::FourBarParameters, q::AbstractVector)
    theta1, theta2 = q
    O = [0.0, 0.0]
    D = [p.l0, 0.0]
    A = p.l1 .* [cos(theta1), sin(theta1)]
    B = A + p.l2 .* [cos(theta2), sin(theta2)]
    G1 = A ./ 2
    G2 = (A + B) ./ 2
    G3 = (D + B) ./ 2
    theta3 = atan(B[2] - D[2], B[1] - D[1])
    return (; O, A, B, D, G1, G2, G3, theta3)
end

"""Phi = (|B-D|^2-l3^2)/2, in m^2."""
function closure_residual(model::FourBarModel, q::AbstractVector)
    f = model.symbolic.functions
    return Float64(f.Phi(q, zeros(2), model.parameter_values)[1])
end

"""
Find a consistent theta2 for a supplied theta1 and assembly branch.
No circle intersections or square roots: solve the scalar *squared-length*
constraint by bisection in an offset angle from the A->D direction.

branch=+1 selects B to the left of the directed line A->D;
branch=-1 selects B to the right. This is only used for initialization.
Afterwards both angles evolve and the constraint is enforced numerically.
"""
function initial_angles(p::FourBarParameters, theta1::Real; branch::Int=1)
    branch in (-1, 1) || throw(ArgumentError("branch must be +1 or -1"))
    A = p.l1 .* [cos(theta1), sin(theta1)]
    dx, dy = p.l0 - A[1], -A[2]     # vector A -> D
    alpha = atan(dy, dx)

    function phi_at(delta)
        theta2 = alpha + branch * delta
        bx = A[1] + p.l2*cos(theta2)
        by = A[2] + p.l2*sin(theta2)
        return ((bx - p.l0)^2 + by^2 - p.l3^2) / 2
    end

    flo, fhi = phi_at(0.0), phi_at(pi)
    scale = max(p.l0^2, p.l1^2, p.l2^2, p.l3^2)
    (flo <= 1e-12*scale && fhi >= -1e-12*scale) ||
        throw(DomainError(theta1, "four-bar cannot close at this input angle"))

    # The squared-length residual increases monotonically on [0, pi].
    lo, hi = 0.0, pi
    for _ in 1:55
        mid = (lo + hi)/2
        if phi_at(mid) > 0
            hi = mid
        else
            lo = mid
        end
    end
    q = [Float64(theta1), alpha + branch*(lo + hi)/2]
    return q
end

"""
Orthogonally project angles onto Phi(q)=0 and angular velocities onto
J(q)*qd=0. Applied at initialization and after accepted ODE steps.
The Newton iteration corrects a scalar constraint, not a circle intersection.
"""
function project_state!(u::AbstractVector, model::FourBarModel;
                        maxiter::Int=12, tol::Float64=1e-12)
    f = model.symbolic.functions
    p = model.parameter_values
    q = @view u[1:2]
    qd = @view u[3:4]
    dummy = zeros(2)
    length_scale = max(model.parameters.l0, model.parameters.l1,
                       model.parameters.l2, model.parameters.l3)
    threshold = tol * length_scale^2

    for _ in 1:maxiter
        phi = Float64(f.Phi(q, dummy, p)[1])
        abs(phi) <= threshold && break
        j = vec(Float64.(f.J(q, dummy, p)))
        jj = dot(j, j)
        jj > 1e-24*length_scale^4 ||
            throw(DomainError(q, "constraint Jacobian is singular"))
        q .-= (phi/jj) .* j
    end
    abs(Float64(f.Phi(q, dummy, p)[1])) <= threshold ||
        throw(ErrorException("holonomic constraint projection failed"))

    j = vec(Float64.(f.J(q, dummy, p)))
    jj = dot(j, j)
    jj > 1e-24*length_scale^4 ||
        throw(DomainError(q, "constraint Jacobian is singular"))
    qd .-= (dot(j, qd)/jj) .* j
    return u
end

"""Reconstruct geometry and output-link angular velocity."""
function configuration(model::FourBarModel, u::AbstractVector)
    q = u[1:2]
    qd = u[3:4]
    pos = body_positions(model.parameters, q)
    p = model.parameters
    vB = p.l1*qd[1] .* [-sin(q[1]), cos(q[1])] +
         p.l2*qd[2] .* [-sin(q[2]), cos(q[2])]
    db = pos.B - pos.D
    omega3 = (db[1]*vB[2] - db[2]*vB[1]) / p.l3^2
    return (; q, qd, positions=pos,
            theta=[q[1], q[2], pos.theta3],
            omega=[qd[1], qd[2], omega3])
end

"""
At each state solve the 3x3 constrained Euler--Lagrange equations:

    [M  J'] [qdd   ] = [Q-c  ]
    [J   0] [lambda]   [-gamma - constraint_stabilization]

The optional Baumgarte term controls drift between projection steps.
"""
function dynamics(model::FourBarModel, u::AbstractVector, t::Real;
                  stabilization_frequency::Real=15.0)
    q, qd = u[1:2], u[3:4]
    f = model.symbolic.functions
    pars = model.parameter_values

    M = Matrix{Float64}(f.M(q, qd, pars))
    c = Vector{Float64}(f.c(q, qd, pars))
    J = Matrix{Float64}(f.J(q, qd, pars))
    gamma = Float64(f.gamma(q, qd, pars)[1])
    phi = Float64(f.Phi(q, qd, pars)[1])

    tau = Float64(model.torque(t))
    Q = [tau, 0.0]

    # Critically damped Baumgarte stabilization of Phi''=0.
    k = Float64(stabilization_frequency)
    k >= 0 || throw(ArgumentError("stabilization_frequency must be nonnegative"))
    constraint_rhs = -gamma - 2*k*dot(vec(J), qd) - k^2*phi

    K = zeros(3, 3)
    K[1:2, 1:2] .= M
    K[1:2, 3] .= vec(J)
    K[3, 1:2] .= vec(J)
    rhs = [Q - c; constraint_rhs]

    x = K \ rhs
    return (; acceleration=x[1:2], lambda=x[3])
end

function rhs!(du, u, model::FourBarModel, t)
    ans = dynamics(model, u, t)
    du[1] = u[3]
    du[2] = u[4]
    du[3] = ans.acceleration[1]
    du[4] = ans.acceleration[2]
    return nothing
end

"""
Simulate [theta1, theta2, omega1, omega2]. With theta2=nothing and
omega2=nothing, derive a constraint-consistent initial state automatically.

This is an acceleration-level (index-reduced) constrained EL solve with
per-accepted-step projection, NOT a direct index-3 DAE integration.
The return value has .t, .u, and .retcode fields.
"""
function simulate(model::FourBarModel;
                  theta1::Real=pi/3, theta2::Union{Nothing,Real}=nothing,
                  omega1::Real=0.0, omega2::Union{Nothing,Real}=nothing,
                  branch::Int=1,
                  tspan::Tuple{<:Real,<:Real}=(0.0, 8.0),
                  fps::Int=30, reltol::Real=1e-9, abstol::Real=1e-11)

    # Use protection, kids.
    fps > 0 || throw(ArgumentError("fps must be positive"))
    tspan[2] > tspan[1] || throw(ArgumentError("tspan must increase"))

    q = isnothing(theta2) ? initial_angles(model.parameters, theta1; branch) :
                            Float64[theta1, theta2]
    f = model.symbolic.functions
    J = vec(Float64.(f.J(q, zeros(2), model.parameter_values)))
    w2 = isnothing(omega2) ?
         (abs(J[2]) > 1e-12 ? -J[1]/J[2]*omega1 :
          throw(DomainError(q, "theta1 is not a valid independent input here"))) :
         omega2
    u0 = Float64[q[1], q[2], omega1, w2]
    project_state!(u0, model)

    problem = ODEProblem(rhs!, u0,
                         (Float64(tspan[1]), Float64(tspan[2])), model)

    # Project the integrator state after each successful adaptive step.
    callback = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> project_state!(integrator.u, model);
        save_positions=(false, false)
    )

    count = max(2, round(Int, (tspan[2]-tspan[1])*fps)+1)
    sampled_times = range(Float64(tspan[1]), Float64(tspan[2]); length=count)
    raw = solve(problem, Tsit5(); callback,
                saveat=sampled_times, save_everystep=false,
                reltol, abstol)

    # saveat can interpolate *inside* an adaptive step. Project those
    # interpolated samples too, so diagnostics/plots are on the manifold.
    states = [project_state!(copy(state), model) for state in raw.u]
    return (; t=Float64.(raw.t), u=states, retcode=raw.retcode)
end

function total_energy(model::FourBarModel, u::AbstractVector)
    return Float64(model.symbolic.functions.energy(
        u[1:2], u[3:4], model.parameter_values))
end

function constraint_multipliers(model::FourBarModel, u::AbstractVector, t::Real)
    return dynamics(model, u, t).lambda
end

end # module FourBarSystem
