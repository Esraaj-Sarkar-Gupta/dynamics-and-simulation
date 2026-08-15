using DifferentialEquations
using Plots

# Define Pendulum
function system!(du, u, p, t)
    # Unpack state vector
    ϕ_1, ϕ_2 = u
    l, g = p

    # Define time differentials
    dϕ_1 = ϕ_2
    dϕ_2 = -g/l * sin(ϕ_1)

    du[1] = dϕ_1
    du[2] = dϕ_2
    return nothing
end

# Setup the problem
u0 = [π/4, 0]
p = [1, 9.810]

tspan = (0.0, 10.0)

# We got a problem buddy?
problem = ODEProblem(system!, u0, tspan, p)
solution = solve(problem, Tsit5())

# State variables plot
p1 = plot(
    solution,
    label = ["Angle" "Angular Velocity"],
    xlabel = "Time (s)",
    ylabel = "State",
)
display(p1)
savefig(p1, "Problem_1/figures/state_plot.png")

# Phase portrait
p2 = plot(
    solution[1, :],
    solution[2, :],
    xlabel = "ϕ₁",
    ylabel = "ϕ₂",
    title = "Phase Space",
    legend = false
)
display(p2)
savefig(p2, "Problem_1/figures/phase_space.png")

# Animate the Pendulum
l = p[1] # We need this

animation = @animate for t in range(tspan[1], tspan[2], length = 200)
    θ = solution(t)[1]
    dθ = solution(t)[2] # We will not be using this result here

    x = l * sin(θ)
    y = -l * cos(θ)

    # Animate the string too!
    plot(
        [0, x],
        [0, y],
        lw = 3,
        color = :black,
        aspect_ratio =:equal,
        xlims=(-1.2l, 1.2l), ylims=(-1.2l, 0.2l)
    )

    # Animate the mass
    scatter!(
        [x], [y], color = :red, ms = 12
    )
end

# It's pronounced gif
gif(animation, "Problem_1/figures/simple_pendulum.gif", fps = 24)