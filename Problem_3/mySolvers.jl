module MySolvers

export ODESolution, Euler, euler_step!

# Define the solution data struct
struct ODESolution
    t::Vector{Float64}
    u::Vector{Vector{Float64}}
end

# Am I the next oiler?
function euler_step!(system!, u, p, t, dt, du)
    system!(du, u, p, t)
    u .+= dt .* du
    return u
end

function Euler(
    system,
    initial_state::Vector{Float64},
    tspan::Tuple{Float64, Float64},
    parameters::Vector{Float64},
    delta_t::Float64
)
    t_start, t_stop = tspan
    number_of_steps = round(Int, (t_stop - t_start) / delta_t)

    t_vec = Vector{Float64}(undef, number_of_steps + 1)
    u_vec = Vector{Vector{Float64}}(undef, number_of_steps + 1)

    u_curr = copy(initial_state)
    du = similar(initial_state)
    t_curr = t_start

    t_vec[1] = t_curr
    u_vec[1] = copy(u_curr)

    for i in 1:number_of_steps
        euler_step!(system, u_curr, parameters, t_curr, delta_t, du)
        t_curr += delta_t

        t_vec[i + 1] = t_curr
        u_vec[i + 1] = copy(u_curr)
    end

    return ODESolution(t_vec, u_vec)
end

end # module