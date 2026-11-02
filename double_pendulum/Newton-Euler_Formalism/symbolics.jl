# ===== Double Pendulum: Symbolic Newton--Euler Derivation ===== #
# Compute the joint moments independently and extract A(θ) and b(θ, ω).
# No numerical solver is required to run this module.

module SymbolicSystem

using Symbolics

export derive_newton_euler, print_symbolic_derivation

# Planar z-component of the cross product.
cross2(a, b) = a[1] * b[2] - a[2] * b[1]

"""
Apply the chain rule to a symbolic expression whose independent variables
are q = [θ1, θ2, ω1, ω2], with q̇ = [ω1, ω2, α1, α2].
"""
function time_derivative(expr, q, qdot)
    terms = [expand_derivatives(Differential(q[i])(expr)) * qdot[i]
             for i in eachindex(q)]
    return simplify(sum(terms))
end

"""
Derive the equations by taking moments about:
  O: the fixed support, using angular momentum of BOTH masses;
  A: the accelerating second joint, using angular momentum of mass 2
     relative to A plus the moving-pivot correction r₂/A × m₂ a_A.

Returns the original moment residuals, reduced residuals, and the
coefficient matrix A and forcing vector b such that A * α = b.
"""
function derive_newton_euler()
    @variables θ1 θ2 ω1 ω2 α1 α2
    @variables m1 m2 l1 l2 g

    q = [θ1, θ2, ω1, ω2]
    qdot = [ω1, ω2, α1, α2]
    dt(expr) = time_derivative(expr, q, qdot)

    # 1. Positions, velocities and accelerations.
    r1 = [l1 * sin(θ1), -l1 * cos(θ1)]
    r2 = [l1 * sin(θ1) + l2 * sin(θ2),
          -l1 * cos(θ1) - l2 * cos(θ2)]

    v1 = dt.(r1)
    v2 = dt.(r2)
    a1 = dt.(v1)

    # 2. Moment about O: dH_O/dt = τ_O.
    HO = cross2(r1, m1 .* v1) + cross2(r2, m2 .* v2)
    τO = -m1 * g * r1[1] - m2 * g * r2[1]
    residual_O = simplify(dt(HO) - τO)

    # 3. Moment about A: dH_A,rel/dt + r₂/A × m₂ a_A = τ_A.
    r2A = r2 .- r1
    v2A = v2 .- v1
    HA_rel = cross2(r2A, m2 .* v2A)
    moving_pivot_term = cross2(r2A, m2 .* a1)
    τA = -m2 * g * r2A[1]
    residual_A = simplify(dt(HA_rel) + moving_pivot_term - τA)

    # 4. Reduce to the two ODEs used in system.jl.
    #    (O - A)/l1 is the θ1 equation; A/(m2*l2) is the θ2 equation.
    eq1 = simplify((residual_O - residual_A) / l1)
    eq2 = simplify(residual_A / (m2 * l2))
    residuals = [eq1, eq2]

    # 5. Because the equations are linear in α1, α2, the Jacobian with
    #    respect to angular accelerations IS the coefficient matrix A.
    accelerations = [α1, α2]
    A = simplify.(Symbolics.jacobian(residuals, accelerations))

    zero_accelerations = Dict(α1 => 0, α2 => 0)
    b = -simplify.([substitute(eq, zero_accelerations) for eq in residuals])

    return (HO=HO, HA_rel=HA_rel,
            moving_pivot_term=moving_pivot_term,
            torque_O=τO, torque_A=τA,
            moment_O=residual_O, moment_A=residual_A,
            reduced_equations=residuals,
            A=A, b=b)
end

"""Print the symbolic derivation and its final matrix formulation."""
function print_symbolic_derivation()
    result = derive_newton_euler()

    println("="^68)
    println("NEWTON--EULER: SYMBOLIC DOUBLE-PENDULUM DERIVATION")
    println("="^68)
    println("\n[1] Angular momentum about fixed pivot O:")
    println("H_O = ", result.HO)
    println("\n[2] Angular momentum relative to moving pivot A:")
    println("H_A,rel = ", result.HA_rel)
    println("Moving-pivot correction r₂/A × m₂ a_A = ", result.moving_pivot_term)
    println("\n[3] Moment equations (each expression = 0):")
    println("O: ", result.moment_O)
    println("A: ", result.moment_A)
    println("\n[4] Reduced equations (each expression = 0):")
    println("Equation 1: ", result.reduced_equations[1])
    println("Equation 2: ", result.reduced_equations[2])
    println("\n[5] Coefficient matrix A(θ) in A * [α1, α2] = b:")
    display(result.A)
    println("\n[6] Forcing vector b(θ, ω):")
    display(result.b)
    println("="^68)

    return result
end

end # module SymbolicSystem

# Allow standalone execution: julia symbolics.jl
if abspath(PROGRAM_FILE) == @__FILE__
    using .SymbolicSystem
    print_symbolic_derivation()
end
