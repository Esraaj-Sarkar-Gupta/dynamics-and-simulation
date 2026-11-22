# Four-bar linkage with Cartesian constraints

Run from the repository root:

```sh
julia four_bar_linkage/fourbar_dae/main.jl
```

This implementation follows the Cartesian DAE derivation in the four-bar
report. The generalized coordinates are
`q = [xA, yA, xB, yB]`, and the state is `[q; q̇]`. The three constraints are
the fixed lengths of `OA`, `AB`, and `DB`.

## Files

- `symbolics.jl` derives the rigid-bar kinetic and potential energies,
  `M_c(q)`, `c_c(q,q̇)`, `J_c(q)`, and `J̇_c(q,q̇)q̇` with Symbolics.jl, then
  compiles numerical functions.
- `system.jl` builds and solves the 7×7 KKT system for `q̈` and the three
  multipliers at every integration step. It also initializes a closed pose
  and tangent velocity, and projects accepted states back onto the position
  and velocity constraints.
- `animations.jl` saves joint trajectories, a snapshot, and an MP4 animation.
- `main.jl` supplies parameters, runs the simulation, and reports constraint
  and energy residuals.

The equations are enforced at acceleration level, so `Tsit5` integrates the
first-order state `(q,q̇)` while the KKT solve determines `(q̈,λ)` at each RHS
evaluation. Position and velocity projection control drift in the original
constraints. This is the index-reduced constrained DAE approach, rather than
passing a residual-form, index-3 system directly to a DAE integrator.

By default, outputs are written to this directory's `output/` folder. Set
`torque=t -> ...` in `main.jl` to apply a driving torque at the input pivot.
