# symbolics.jl
# Symbolic derivation of the Mass Matrix M(θ) and Forcing Vector f(θ, dθ)

using Symbolics

function derive_triple_pendulum_equations()
    println("1. Setting up symbolic framework...")
    @variables t
    d = Differential(t)
    dd = Differential(t)^2

    @variables m1 m2 m3 l1 l2 l3 g
    @variables θ1(t) θ2(t) θ3(t)
    @variables q1 q2 q3 dq1 dq2 dq3

    # Kinematics and Energy
    V = -g * (m1*l1*cos(q1) + m2*(l1*cos(q1) + l2*cos(q2)) + m3*(l1*cos(q1) + l2*cos(q2) + l3*cos(q3)))

    # Unsigned velocities (using simple chain-rule differentiation)
    vx1 = l1 * cos(q1) * dq1
    vy1 = l1 * sin(q1) * dq1

    vx2 = vx1 + l2 * cos(q2) * dq2
    vy2 = vy1 + l2 * sin(q2) * dq2

    vx3 = vx2 + l3 * cos(q3) * dq3
    vy3 = vy2 + l3 * sin(q3) * dq3

    T = 0.5 * m1 * (vx1^2 + vy1^2) + 0.5 * m2 * (vx2^2 + vy2^2) + 0.5 * m3 * (vx3^2 + vy3^2)
    L = T - V

    # Partials evaluated using static variables
    sub_back = Dict(q1=>θ1, q2=>θ2, q3=>θ3, dq1=>d(θ1), dq2=>d(θ2), dq3=>d(θ3))

    pL_pθ1 = substitute(expand_derivatives(Differential(q1)(L)), sub_back)
    pL_pθ2 = substitute(expand_derivatives(Differential(q2)(L)), sub_back)
    pL_pθ3 = substitute(expand_derivatives(Differential(q3)(L)), sub_back)

    pL_pdθ1 = substitute(expand_derivatives(Differential(dq1)(L)), sub_back)
    pL_pdθ2 = substitute(expand_derivatives(Differential(dq2)(L)), sub_back)
    pL_pdθ3 = substitute(expand_derivatives(Differential(dq3)(L)), sub_back)

    # Euler-Lagrange equations (LHS = 0)
    eq1_lhs = expand_derivatives(d(pL_pdθ1)) - pL_pθ1
    eq2_lhs = expand_derivatives(d(pL_pdθ2)) - pL_pθ2
    eq3_lhs = expand_derivatives(d(pL_pdθ3)) - pL_pθ3

    println("2. Extracting Mass matrix and forcing vectors...")
    vars_dd = [dd(θ1), dd(θ2), dd(θ3)]
    eqs_lhs = [eq1_lhs, eq2_lhs, eq3_lhs]

    M_sym = Symbolics.jacobian(eqs_lhs, vars_dd)

    zero_dd = Dict(dd(θ1)=>0, dd(θ2)=>0, dd(θ3)=>0)
    f_sym = substitute.(-1 .* eqs_lhs, (zero_dd,))

    println("3. Compiling to fast numerical functions...")
    @variables th1 th2 th3 w1 w2 w3
    state_subs = Dict(θ1=>th1, θ2=>th2, θ3=>th3, d(θ1)=>w1, d(θ2)=>w2, d(θ3)=>w3)

    M_clean = substitute.(M_sym, (state_subs,))
    f_clean = substitute.(f_sym, (state_subs,))

    state_vars = [th1, th2, th3, w1, w2, w3]
    param_vars = [m1, m2, m3, l1, l2, l3, g]

    M_func = build_function(M_clean, state_vars, param_vars, expression=Val{false})[1]
    f_func = build_function(f_clean, state_vars, param_vars, expression=Val{false})[1]

    return M_func, f_func
end