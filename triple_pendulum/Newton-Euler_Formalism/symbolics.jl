# ===== Triple Pendulum: Symbolic Newton--Euler Derivation ===== #
# Derive the moment equations at the fixed joint O, moving joint A (at m1),
# and moving joint B (at m2). Extract a 3x3 matrix A and a forcing vector b.

module SymbolicSystem

using Symbolics

export derive_newton_euler, print_symbolic_derivation, validate_symbolic_matrix

# Signed z-component of the cross product of two planar vectors.
cross2(a, b) = a[1]*b[2] - a[2]*b[1]

"""Differentiate a function of [θ1,θ2,θ3,ω1,ω2,ω3] in time by the chain rule."""
function time_derivative(expr, state, state_dot)
    # The fixed support is represented by integer zeros.
    expr isa Integer && return zero(expr)
    terms = [expand_derivatives(Differential(state[i])(expr))*state_dot[i]
             for i in eachindex(state)]
    return sum(terms)
end

"""
Derive the Newton--Euler dynamics without using the Lagrangian.

For pivot P_{j-1}, take moments of the subchain of bobs j,...,3:
    E_j = dH_{j,rel}/dt + Σ_{k=j}^3 (r_k-r_P) × m_k*a_P - τ_{g,j} = 0.

The moving-pivot term is mandatory when j=2,3. Define reduced equations
E1-E2=0, E2-E3=0, E3=0, then extract the acceleration coefficients using
Symbolics.jacobian. The result is A(θ)*α=b(θ,ω).
"""
function derive_newton_euler(; simplify_output=false)
    @variables θ1 θ2 θ3 ω1 ω2 ω3 α1 α2 α3
    @variables m1 m2 m3 l1 l2 l3 g

    θ = [θ1, θ2, θ3]
    ω = [ω1, ω2, ω3]
    α = [α1, α2, α3]
    m = [m1, m2, m3]
    l = [l1, l2, l3]
    state = vcat(θ, ω)
    state_dot = vcat(ω, α)
    dt(expr) = time_derivative(expr, state, state_dot)

    # r[1],v[1],a[1] are the inertially fixed support O.
    # r[2:4] correspond to masses m1,m2,m3.
    r = Any[[0, 0]]
    for j in 1:3
        previous = r[end]
        push!(r, [previous[1]+l[j]*sin(θ[j]),
                  previous[2]-l[j]*cos(θ[j])])
    end

    # Differentiating positions *symbolically* with the full chain rule.
    v = [dt.(rj) for rj in r]
    a = [dt.(vj) for vj in v]

    H_relative = []
    dH_relative = []
    correction = []
    torque_gravity = []
    residuals = []

    for j in 1:3
        # j=1: whole chain about O; j=2: (m2,m3) about A;
        # j=3: m3 about B.
        rp, vp, ap = r[j], v[j], a[j]
        H_j = 0
        dH_j = 0
        C_j = 0
        τ_j = 0
        for k in j:3
            rk = r[k+1]
            vk = v[k+1]
            relr = rk .- rp
            relv = vk .- vp
            rela = a[k+1] .- ap
            H_j += m[k]*cross2(relr, relv)
            # d/dt(relr × m*relv) = relr × m*rela because
            # relv × relv = 0. This avoids unnecessary product expansion.
            dH_j += m[k]*cross2(relr, rela)
            C_j += m[k]*cross2(relr, ap)
            τ_j -= m[k]*g*relr[1]
        end
        # No joint reaction torque: the force at the chosen pivot is radial.
        push!(H_relative, H_j)
        push!(dH_relative, dH_j)
        push!(correction, C_j)
        push!(torque_gravity, τ_j)
        push!(residuals, dH_j+C_j-τ_j)
    end

    # E1-E2 isolates θ1, E2-E3 isolates θ2, E3 isolates θ3.
    equations = [residuals[1]-residuals[2],
                 residuals[2]-residuals[3],
                 residuals[3]]

    # Acceleration linearity -> the Jacobian IS the mass matrix.
    M_sym = Symbolics.jacobian(equations, α)
    zero_alpha = Dict(α1=>0, α2=>0, α3=>0)
    b_sym = [-substitute(E, zero_alpha) for E in equations]
    # Full trig simplification is optional: it can be costly for three links.
    if simplify_output
        M_sym = simplify.(M_sym)
        b_sym = simplify.(b_sym)
    end

    return (positions=r, velocities=v, accelerations=a,
            angular_momenta=H_relative, angular_momentum_derivatives=dH_relative,
            moving_pivot_terms=correction,
            gravity_torques=torque_gravity, raw_moment_equations=residuals,
            reduced_equations=equations, M=M_sym, b=b_sym,
            angles=θ, omegas=ω, alphas=α, masses=m, lengths=l, gravity=g)
end

"""
Compare symbolic M and b numerically with the compact Newton--Euler
coefficients implemented in system.jl (a cheap, independent spot check).
No reference to System is needed; expected values are written explicitly.
"""
function validate_symbolic_matrix(result; atol=1e-8)
    θ = [0.31, -0.47, 0.86]
    ω = [0.20, -0.36, 0.55]
    m = [1.2, 0.8, 1.4]
    l = [0.7, 1.1, 0.9]
    g = 9.81
    subs = Dict{Any,Any}()
    for (key, val) in zip(vcat(result.angles, result.omegas,
                               result.masses, result.lengths, [result.gravity]),
                           vcat(θ, ω, m, l, [g]))
        subs[key] = val
    end

    # Numeric substitution must fold arithmetic before unwrapping Num.
    # Float64(::Symbolics.Num) fails when the symbolic wrapper remains.
    function evaluate_numeric(expr)
        evaluated = substitute(expr, subs; fold=Val(true))
        numeric = evaluated isa Symbolics.Num ? Symbolics.value(evaluated) : evaluated
        numeric = numeric isa Number && !(numeric isa Symbolics.Num) ? numeric :
                  Symbolics.symbolic_to_float(numeric)
        numeric isa Real && !(numeric isa Symbolics.Num) ||
            error("Unresolved symbols after substitution: $(evaluated)")
        return Float64(numeric)
    end

    M_num = map(evaluate_numeric, result.M)
    b_num = map(evaluate_numeric, result.b)
    S = [sum(m[i:3]) for i in 1:3]
    M_expected = [S[max(i,j)]*l[i]*l[j]*cos(θ[i]-θ[j])
                  for i in 1:3, j in 1:3]
    b_expected = [-S[i]*g*l[i]*sin(θ[i]) -
                  sum(S[max(i,j)]*l[i]*l[j]*sin(θ[i]-θ[j])*ω[j]^2
                      for j in 1:3 if j != i)
                  for i in 1:3]
    matrix_error = maximum(abs.(M_num .- M_expected))
    rhs_error = maximum(abs.(b_num .- b_expected))
    println("Symbolic-vs-compact max |ΔM| = ", matrix_error)
    println("Symbolic-vs-compact max |Δb| = ", rhs_error)
    (matrix_error < atol && rhs_error < atol) ||
        error("The symbolic joint-moment equations disagree with the numerical model.")
    return true
end

"""Print derivation ingredients and the extracted 3x3 Newton--Euler system."""
function print_symbolic_derivation(; validate=true, verbose=false)
    result = derive_newton_euler() # no expensive global simplification
    println("="^72)
    println("TRIPLE PENDULUM: NEWTON--EULER MOMENTS AT O, A, B")
    println("="^72)
    for (j, label) in enumerate(("O (all three masses)",
                                 "A (masses 2 and 3)",
                                 "B (mass 3)"))
        println("\n[$j] Moment about $label:")
        if verbose
            println("   H_rel = ", result.angular_momenta[j])
            println("   moving-pivot correction = ", result.moving_pivot_terms[j])
            println("   gravitational torque = ", result.gravity_torques[j])
        else
            println("   Angular momentum, pivot correction and gravity torque derived.")
        end
        # Large raw equations can overwhelm the terminal; inspect result.raw_moment_equations.
    end
    println("\nReduced residuals: [E_O-E_A, E_A-E_B, E_B] = 0")
    println("\nMass matrix M(θ) extracted with jacobian(E_reduced, α):")
    if verbose
        display(result.M)
    else
        println("   3 × 3 symbolic matrix constructed (use verbose=true to display).")
    end
    println("\nForcing vector b(θ,ω) = -E_reduced|_{α=0}:")
    if verbose
        display(result.b)
    else
        println("   3-element symbolic vector constructed (use verbose=true to display).")
    end
    if validate
        println("\nIndependent numerical consistency check:")
        validate_symbolic_matrix(result)
    end
    return result
end

end # module SymbolicSystem

if abspath(PROGRAM_FILE) == @__FILE__
    using .SymbolicSystem
    print_symbolic_derivation()
end
