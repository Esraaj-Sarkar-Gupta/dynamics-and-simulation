# Define a Lorentz System
function Lorentz_system!(du, u, p, t)
    # Unpack the system
    x, y, z = u
    σ, ρ, β = p

    # Dynamics
    dx = σ*(y - x)
    dy = x*(ρ - z) - y
    dz = x*y - β*z

    du[1] = dx; du[2] = dy; du[3] = dz
    return nothing # Absolutely nothing at all!
end