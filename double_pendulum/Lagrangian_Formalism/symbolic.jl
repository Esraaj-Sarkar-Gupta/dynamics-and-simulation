# ===== Dynamics and Simulations ===== #
# October 2026
# Esraaj Sarkar Gupta

# Symbolic simplification of the double pendulum EOM
# from the Lagrangian.

using Symbolics

# Define variables
@variables t
D = Differential(t)

@variables m1 m2 l1 l2 g

# Define functions of time
@variables θ1(t) θ2(t)

# Define derivatives
θ1_dot = D(θ1)
θ2_dot = D(θ2)

θ1_ddot = D(θ1_dot)
θ2_ddot = D(θ2_dot)

# Define static variables to evaluate partial derivatives
@variables q1 q2 q1_dot q2_dot

# Define static varibles to eval partial derivatives
L_static = 0.5 * (m1 + m2) * l1^2 * q1_dot^2 + 
           0.5 * m2 * l2^2 * q2_dot^2 + 
           m2 * l1 * l2 * q1_dot * q2_dot * cos(q1 - q2) + 
           (m1 + m2) * g * l1 * cos(q1) + 
           m2 * g * l2 * cos(q2)

# Evaulate partial differentials
dL_dq1_dot = expand_derivatives(
    Differential(q1_dot)(L_static)
)

dL_dq2_dot = expand_derivatives(
    Differential(q2_dot)(L_static)
)

dL_dq1 = expand_derivatives(
    Differential(q1)(L_static)
)

dL_dq2 = expand_derivatives(
    Differential(q2)(L_static)
)

# Substitute the time-dependent dynamic functions into the expressions
sub_back = Dict(
    q1 => θ1,
    q2 => θ2,
    q1_dot => θ1_dot,
    q2_dot => θ2_dot
)

pL_pθ1_dot = substitute(dL_dq1_dot, sub_back)
pL_pθ2_dot = substitute(dL_dq2_dot, sub_back)
pL_pθ1 = substitute(dL_dq1, sub_back)
pL_pθ2 = substitute(dL_dq2, sub_back)

# Form the Euler-Lagrange Equation
eq1 = expand_derivatives(D(pL_pθ1_dot)) - pL_pθ1 ~ 0
eq2 = expand_derivatives(D(pL_pθ2_dot)) - pL_pθ2 ~ 0

# Solve
solutions = Symbolics.solve_for(
    [eq1, eq2], [θ1_ddot, θ2_ddot]
)

# Simplify the final explicit equations of motion
θ1_ddot_explicit = simplify(solutions[1])
θ2_ddot_explicit = simplify(solutions[2])

println("Equation for θ1_ddot:")
println(θ1_ddot_explicit)
println("\nEquation for θ2_ddot:")
println(θ2_ddot_explicit)