module FourBarSymbolics

using Symbolics

export derive_euler_lagrange

"""
Symbolically derive the constrained two-angle Euler--Lagrange equations:

    M(q) * qdd + c(q, qd) + J(q)' * lambda = [tau, 0]
    Phi(q) = 0

q = [theta1, theta2] are the absolute angles of OA and AB.
B = A + l2*[cos(theta2), sin(theta2)] is constructed directly.
The single constraint is Phi = (|B-D|^2-l3^2)/2.

The third bar's angular velocity is reconstructed from B's velocity,
using the fact that |B-D|=l3 and (B-D) dot Bdot=0 on the constraint.

All derivatives are produced by Symbolics.jl. No pre-expanded
Euler--Lagrange equations, hand-derived mass matrices, or simplify calls.
Numerical functions take (q, qd, parameters), with parameters ordered
[l0, l1, l2, l3, m1, m2, m3, g].
"""
function derive_euler_lagrange()
    @variables theta1 theta2 omega1 omega2
    @variables l0 l1 l2 l3 m1 m2 m3 g

    q = [theta1, theta2]
    qd = [omega1, omega2]
    pars = [l0, l1, l2, l3, m1, m2, m3, g]

    # Construct O -> A -> B directly; D is the other fixed pivot.
    O = [0, 0]
    D = [l0, 0]

    e1 = [cos(theta1), sin(theta1)]
    e2 = [cos(theta2), sin(theta2)]
    n1 = [-sin(theta1), cos(theta1)]
    n2 = [-sin(theta2), cos(theta2)]

    A = l1 .* e1
    B = A + l2 .* e2
    DB = B - D

    vA = l1 * omega1 .* n1
    vB = vA + l2 * omega2 .* n2

    # Centres of mass of uniform bars (including the ground-pinned DB bar).
    G1 = A ./ 2
    G2 = (A + B) ./ 2
    G3 = (D + B) ./ 2
    vG1 = vA ./ 2
    vG2 = (vA + vB) ./ 2
    vG3 = vB ./ 2

    cross2(a, b) = a[1]*b[2] - a[2]*b[1]
    norm2(a) = sum(x -> x^2, a)

    I1 = m1*l1^2/12
    I2 = m2*l2^2/12
    I3 = m3*l3^2/12

    # On the constraint manifold: omega3 = (DB x vB)/l3^2.
    omega3 = cross2(DB, vB) / l3^2

    T = (m1/2)*norm2(vG1) + (I1/2)*omega1^2 +
        (m2/2)*norm2(vG2) + (I2/2)*omega2^2 +
        (m3/2)*norm2(vG3) + (I3/2)*omega3^2

    V = g * (m1*G1[2] + m2*G2[2] + m3*G3[2])
    L = T - V

    Phi = [(norm2(DB) - l3^2) / 2]

    # Terms of the Lagrangian, derived from the chain rule
    momentum = [Symbolics.derivative(L, v; simplify=false) for v in qd]
    M = Symbolics.jacobian(momentum, qd; simplify=false)

    dmomentum_dq = Symbolics.jacobian(momentum, q; simplify=false)
    dL_dq = [Symbolics.derivative(L, qi; simplify=false) for qi in q]

    c = dmomentum_dq * qd - dL_dq

    J = Symbolics.jacobian(Phi, q; simplify=false)  # 1 x 2
    gamma = [sum(Symbolics.derivative(J[1, i], q[j]; simplify=false) *
                 qd[i] * qd[j] for i in 1:2 for j in 1:2)]

    # Arrays return (out-of-place, in-place); use the first one.
    compile_array(x) = Symbolics.build_function(
        x, q, qd, pars; expression=Val{false}
    )[1]
    # Scalars return the out-of-place function directly.
    compile_scalar(x) = Symbolics.build_function(
        x, q, qd, pars; expression=Val{false}
    )

    functions = (
        M      = compile_array(M),
        c      = compile_array(c),
        J      = compile_array(J),
        gamma  = compile_array(gamma),
        Phi    = compile_array(Phi),
        energy = compile_scalar(T + V),
    )

    expressions = (; A, B, DB, G1, G2, G3, omega3,
                    T, V, L, Phi, M, c, J, gamma)
    return (; expressions, functions)
end

end # module FourBarSymbolics
