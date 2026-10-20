import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time as tm

def mag(x): # Define magnitude function
    return abs(x)

# Constants
rho = 1.74 # This value varies from mass 1 and mass 2, but we'll assume it doesn't
l1 = 0.245
l2 = 0.1778
m1 = 0.4263
m2 = 0.15486
g = 9.81

# The system of ODEs
def equations(t, y):
    theta1, omega1, theta2, omega2 = y
    
    # Coefficients
    A11 = (1/3) * rho * l1**3 + rho * l1**2 * l2
    A12 = (1/2) * l1 * l2**2 * rho
    A21 = (1/2) * l1 * l2**2 * rho
    A22 = (1/3) * l2**3 * rho

    # Right-hand side (gravitational forces)
    b1 = -g * l1 * np.sin(theta1) * (m1/2 + m2)
    b2 = -g * m2 * l2 * np.sin(theta2) / 2
    
    # Solve for angular accelerations (ddot_theta1, ddot_theta2)
    A = np.array([[A11, A12], [A21, A22]])
    b = np.array([b1, b2])
    
    ddot_theta = np.linalg.solve(A, b)
    
    return [omega1, ddot_theta[0], omega2, ddot_theta[1]]

# Initial conditions: [theta1, omega1, theta2, omega2]
y0 = [np.pi/2, 0, 0, 0]  # initial angles and angular velocities

# Time span
t_span = (0, 50)
t_eval = np.linspace(*t_span, 1000)

# Solve the system of ODEs
start_time = tm.time()
solution = solve_ivp(equations, t_span, y0, t_eval=t_eval)
end_time = tm.time()
print(f"> System of equations solved in {end_time - start_time:.4f} seconds.")

# Plot the results
plt.figure(figsize=(10, 5))

plt.plot(solution.t, solution.y[0], label='theta1 (rad)', color='green')
plt.plot(solution.t, solution.y[2], label='theta2 (rad)', color='red')

plt.title('Double Pendulum Angles over Time')
plt.xlabel('Time (s)')
plt.ylabel('Angle (rad)')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 5))

plt.plot(solution.y[0], solution.y[1], label='Mass 1', color='green')
plt.plot(solution.y[2], solution.y[3], label='Mass 2', color='red')

plt.title("Double Pendulum Phase Space")
plt.xlabel('Angle (rad)')
plt.ylabel('Angular velocity (rad/s)')
plt.legend()
plt.grid(True)
plt.show()

# Setup Animation
start_time_2 = tm.time()

fig, ax = plt.subplots(figsize=(10, 10))
d_limit = mag(l1) + mag(l2) + 0.5
ax.set_xlim(-d_limit, d_limit)
ax.set_ylim(-d_limit, d_limit)
ax.set_aspect('equal')
ax.grid(True)

# Static origin marker
ax.plot(0, 0, 'ko', markersize=8)

# Dynamic line and scatter elements
line1, = ax.plot([], [], linewidth=5, color='green')
line2, = ax.plot([], [], linewidth=5, color='red')
cog1, = ax.plot([], [], 'o', color='blue')
cog2, = ax.plot([], [], 'o', color='blue')
trace, = ax.plot([], [], color='purple')

cog2_x, cog2_y = [], []

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    cog1.set_data([], [])
    cog2.set_data([], [])
    trace.set_data([], [])
    return line1, line2, cog1, cog2, trace

def update(frame):
    theta1 = solution.y[0][frame]
    theta2 = solution.y[2][frame]

    # Positions of rod endpoints
    x1 = l1 * np.sin(theta1)
    y1 = -l1 * np.cos(theta1)

    x2 = x1 + l2 * np.sin(theta2)
    y2 = y1 - l2 * np.cos(theta2)

    # Centers of gravity
    c1_x = 0.5 * l1 * np.sin(theta1)
    c1_y = -0.5 * l1 * np.cos(theta1)

    c2_x = x1 + 0.5 * l2 * np.sin(theta2)
    c2_y = y1 - 0.5 * l2 * np.cos(theta2)

    cog2_x.append(c2_x)
    cog2_y.append(c2_y)

    line1.set_data([0, x1], [0, y1])
    line2.set_data([x1, x2], [y1, y2])
    cog1.set_data([c1_x], [c1_y])
    cog2.set_data([c2_x], [c2_y])
    trace.set_data(cog2_x, cog2_y)

    return line1, line2, cog1, cog2, trace

ani = animation.FuncAnimation(
    fig, update, frames=len(solution.t), init_func=init, blit=True, interval=33
)

print("> Rendering and saving 'dpencu_c.gif'...")
ani.save('dpencu_c.gif', writer='pillow', fps=30)
plt.close(fig)

end_time_2 = tm.time()
print(f"> Finished saving 'dpencu_c.gif' in current directory in {end_time_2 - start_time_2:.4f} seconds.")

# Plot final path covered by second mass
plt.figure(figsize=(10, 10))
plt.scatter(0, 0, color='red', s=50)
plt.plot(cog2_x, cog2_y, color='purple')
plt.title("System Path - Path covered by second mass")
plt.grid(True)
plt.show()