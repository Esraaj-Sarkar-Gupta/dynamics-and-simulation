module System

export DoublePenduParams, double_pendulum!

Base.@kwdef struct DoublePenduParams
    m1 :: Float64 = 1.0
    m2 :: Float64 = 1.0
    l1 :: Float64 = 1.0
    l2 :: Float64 = 1.0
    g  :: Float64 = 9.810
end

function double_pendulum!(du, u, p::DoublePenduParams, t)
    # Unpack state: positions and velocities
    x1, y1, x2, y2, vx1, vy1, vx2, vy2 = u
    
    # Relative positions and velocities
    Δx = x2 - x1
    Δy = y2 - y1
    Δvx = vx2 - vx1
    Δvy = vy2 - vy1
    
    m1, m2, g = p.m1, p.m2, p.g

    # Initialize 10x10 matrix A and 10x1 vector b
    A = zeros(eltype(u), 10, 10)
    b = zeros(eltype(u), 10)

    # 1-4: Kinematic mappings
    A[1, 1] = 1.0
    A[2, 2] = 1.0
    A[3, 3] = 1.0
    A[4, 4] = 1.0
    
    b[1] = vx1
    b[2] = vy1
    b[3] = vx2
    b[4] = vy2

    # 5-8: Newton's Second Law
    A[5, 5] = m1
    A[5, 9] = 2 * x1
    A[5, 10] = -2 * Δx
    b[5] = 0.0

    A[6, 6] = m1
    A[6, 9] = 2 * y1
    A[6, 10] = -2 * Δy
    b[6] = -m1 * g

    A[7, 7] = m2
    A[7, 10] = 2 * Δx
    b[7] = 0.0

    A[8, 8] = m2
    A[8, 10] = 2 * Δy
    b[8] = -m2 * g

    # 9-10: Acceleration-level constraints
    A[9, 5] = x1
    A[9, 6] = y1
    b[9] = -(vx1^2 + vy1^2)

    A[10, 5] = -Δx
    A[10, 6] = -Δy
    A[10, 7] = Δx
    A[10, 8] = Δy
    b[10] = -(Δvx^2 + Δvy^2)

    # Solve the linear system
    z = A \ b

    # Update state derivatives in-place
    du[1:8] .= z[1:8]
    
    return nothing
end

end # module