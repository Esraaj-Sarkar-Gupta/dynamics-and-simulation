module FourBarDAESymbolics

using Symbolics

export derive_cartesian_dae

"""Derive Cartesian constrained dynamics for the three moving uniform bars.

Coordinate order is q = [xA, yA, xB, yB], with velocity order matching q.
Parameter order is [l0, l1, l2, l3, m1, m2, m3, g].
"""
function derive_cartesian_dae()
    @variables xA yA xB yB vxA vyA vxB vyB
    @variables l0 l1 l2 l3 m1 m2 m3 g

    q = [xA, yA, xB, yB]
    qd = [vxA, vyA, vxB, vyB]
    pars = [l0, l1, l2, l3, m1, m2, m3, g]

    rA, rB, rD = [xA, yA], [xB, yB], [l0, 0]
    vA, vB = [vxA, vyA], [vxB, vyB]
    dAB, dDB = rB - rA, rB - rD
    vAB = vB - vA

    cross2(a, b) = a[1]*b[2] - a[2]*b[1]
    norm2(a) = sum(x -> x^2, a)

    I1, I2, I3 = m1*l1^2/12, m2*l2^2/12, m3*l3^2/12
    omega1 = cross2(rA, vA) / l1^2
    omega2 = cross2(dAB, vAB) / l2^2
    omega3 = cross2(dDB, vB) / l3^2

    # Midpoint velocities of each uniform moving link.
    vG1 = vA / 2
    vG2 = (vA + vB) / 2
    vG3 = vB / 2
    T = m1*norm2(vG1)/2 + m2*norm2(vG2)/2 + m3*norm2(vG3)/2 +
        (I1*omega1^2 + I2*omega2^2 + I3*omega3^2)/2
    V = (m1*g*yA + m2*g*(yA + yB) + m3*g*yB)/2
    L = T - V

    Psi = [
        (norm2(rA) - l1^2)/2,
        (norm2(dAB) - l2^2)/2,
        (norm2(dDB) - l3^2)/2,
    ]
    J = Symbolics.jacobian(Psi, q; simplify=false)

    momentum = [Symbolics.derivative(L, v; simplify=false) for v in qd]
    M = Symbolics.jacobian(momentum, qd; simplify=false)
    c = Symbolics.jacobian(momentum, q; simplify=false) * qd -
        [Symbolics.derivative(L, qi; simplify=false) for qi in q]
    Jdot = [sum(Symbolics.derivative(J[a, i], q[j]; simplify=false) * qd[j]
                 for j in 1:4) for a in 1:3, i in 1:4]
    gamma = Jdot * qd

    compile_array(x) = Symbolics.build_function(
        x, q, qd, pars; expression=Val{false}
    )[1]
    compile_scalar(x) = Symbolics.build_function(
        x, q, qd, pars; expression=Val{false}
    )

    functions = (
        M = compile_array(M),
        c = compile_array(c),
        J = compile_array(J),
        Jdot = compile_array(Jdot),
        gamma = compile_array(gamma),
        Psi = compile_array(Psi),
        energy = compile_scalar(T + V),
    )
    expressions = (; q, qd, pars, T, V, L, Psi, M, c, J, Jdot, gamma,
                    omega1, omega2, omega3)
    return (; expressions, functions)
end

end # module FourBarDAESymbolics
