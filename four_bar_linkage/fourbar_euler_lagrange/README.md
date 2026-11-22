# Four-bar linkage — two-angle constrained Euler–Lagrange in Julia

## Layout

- `main.jl`: one entry point; parameters, initialization, diagnostics, plots.
- `symbolics.jl`: derive `T`, `V`, `L`, `M`, `c`, `J`, `gamma` with `Symbolics.jl` and compile callable functions; no hand-expanded EL equations or `simplify`.
- `system.jl`: directly construct A and B from absolute angles `theta1`, `theta2`, enforce the DB constraint, integrate dynamics.
- `animations.jl`: save `output/trajectories.png`, `output/snapshot.png`, and `output/fourbar_euler_lagrange.mp4`.

## Install

With Julia installed (Julia 1.10 or later recommended), in this directory:

```bash
julia --project=. -e 'using Pkg; Pkg.add(["Symbolics", "DifferentialEquations", "CairoMakie"])'
julia --project=. main.jl
```

The first command creates `Project.toml` and `Manifest.toml` in this directory.

## Geometry

`O=(0,0)`, `D=(l0,0)` remain fixed. The forward kinematics need no square roots:

```text
A = l1 [cos(theta1), sin(theta1)]
B = A + l2 [cos(theta2), sin(theta2)]
DB = B - D
theta3 = atan(DB_y, DB_x)
```

`theta2` is *not* freely specifiable: the final bar must retain its length:

```text
Phi(theta1,theta2) = (|B-D|^2-l3^2)/2 = 0
```

There is one holonomic constraint for two angles, leaving one degree of freedom.

For a uniform bar, each mass is distributed (no point-mass approximation).
The output link's angular speed is `omega3 = (DB cross vB)/l3^2`, valid
on the constraint manifold.

## Dynamics

`Symbolics.jl` differentiates the Lagrangian to obtain

```text
M(q) qdd + c(q,qd) + J(q)' lambda = [tau(t), 0]
J(q) qdd = -gamma(q,qd)
```

The integrator solves the 3x3 block system for both angular accelerations
and the constraint multiplier. It integrates
`u=[theta1,theta2,omega1,omega2]` with `Tsit5` from DifferentialEquations.
It uses a damped acceleration-level constraint to reduce numerical drift,
and orthogonal projection after every accepted step to restore `Phi=0` and
`J*qd=0`. The interpolated output frames are separately projected.
This is an **index-reduced constrained EL simulation with projection**,
not a direct index-3 DAE integrator.

An admissible theta2 for the supplied initial theta1 is found *once*
by bisection of the squared-length constraint, selecting the `branch`.
This is not a circle-intersection calculation and uses no square roots.
Branch choice affects initialization only; the two angles evolve together.

## Outputs and checks

- `output/trajectories.png`: paths A/B, angles 1/2/3, angular speeds,
  and total mechanical energy.
- `output/snapshot.png`: linkage at the middle saved time.
- `output/fourbar_euler_lagrange.mp4`: the animated linkage.

For zero torque, inspect mechanical energy drift and `|DB|-l3`.
Geometrically singular configurations (zero constraint Jacobian) cannot
be handled by this coordinate chart. A finite angular input torque changes
mechanical energy, so energy conservation then does not apply.

## Validation

The mathematical implementation and API calls have been reviewed, but the
Julia files have not been executed in this environment (Julia unavailable).
Run `main.jl` in your Julia environment to verify exact package compatibility.
