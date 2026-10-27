# symbolics.jl
# Symbolic derivation of the Mass Matrix M(θ) and Forcing Vector f(θ, dθ)

using Symbolics

# 1. Define independent variable (time) and differential operator
@variables t
d = Differential(t)

# 2. Define time-dependent position variables and constants
@variables x1(t) y1(t) x2(t) y2(t) x3(t) y3(t)
@variables l1 l2 l3

# Relative displacements
Δx2 = x2 - x1
Δx3 = x3 - x2
Δy2 = y2 - y1
Δy3 = y3 - y2

# Position constraints
Φ1 = x1^2 + y1^2 - l1^2
Φ2 = Δx2^2 + Δy2^2 - l2^2
Φ3 = Δx3^2 + Δy3^2 - l3^2

# 3. Compute second time derivatives (expand_derivatives executes the chain rule)
ddΦ1 = simplify(expand_derivatives(d(d(Φ1))))
ddΦ2 = simplify(expand_derivatives(d(d(Φ2))))
ddΦ3 = simplify(expand_derivatives(d(d(Φ3))))

# Print position constraints
println("="^50)
println("POSITION CONSTRAINTS [ Φ(q) = 0 ]")
println("="^50)
println("Φ1 = ", Φ1)
println("Φ2 = ", Φ2)
println("Φ3 = ", Φ3)

# Print second time derivatives (acceleration constraints)
println("\n" * "="^50)
println("ACCELERATION CONSTRAINTS [ d²Φ/dt² = 0 ]")
println("="^50)
println("d²Φ1/dt² = ", ddΦ1)
println("\nd²Φ2/dt² = ", ddΦ2)
println("\nd²Φ3/dt² = ", ddΦ3)