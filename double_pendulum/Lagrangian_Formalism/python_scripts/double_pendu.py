print("Simulating a double pendulum...")

import numpy as np
from scipy.integrate import odeint
import time as tm
import matplotlib.pyplot as plt
import matplotlib.animation as animation

print("- Imports complete")

def sys(state, t, m1, m2, l1, l2, g):
    theta1, theta2, theta1_dot, theta2_dot = state

    theta1_ddot = (
        -g * (2*m1 + m2) * np.sin(theta1)
        - m2 * g * np.sin(theta1 - 2*theta2)
        - 2 * np.sin(theta1 - theta2) * m2 * (
            l2 * (theta2_dot**2) + l1 * (theta1_dot**2) * np.cos(theta1 - theta2)
        )
    ) / (l1 * (2*m1 + m2 - m2 * np.cos(2*theta1 - 2*theta2)))

    theta2_ddot = (
        2 * np.sin(theta1 - theta2) * (
            l1 * (theta1_dot**2) * (m1 + m2)
            + g * (m1 + m2) * np.cos(theta1)
            + (theta2_dot**2) * l2 * m2 * np.cos(theta1 - theta2)
        )
    ) / (l2 * (2*m1 + m2 - m2 * np.cos(2*theta1 - 2*theta2)))
    return np.array([theta1_dot, theta2_dot, theta1_ddot, theta2_ddot])

print("- System function defined")

# Initial states:
theta1_0 = 0
theta1_dot_0 = 1
theta2_0 = np.pi / 2
theta2_dot_0 = 1

state0 = np.array([theta1_0, theta2_0, theta1_dot_0, theta2_dot_0])
print("> Initial states defined")

# System parameters:
m1 = 10
m2 = 5
l1 = 4
l2 = 10
g = 9.81 * 0

T = 20 # Max time (seconds)
fps = 30
P = T * fps # Number of data points

t = np.linspace(0, T, P)
print("> System parameters defined")

s_time = tm.time()
sol = odeint(sys, state0, t, args=(m1, m2, l1, l2, g))
e_time = tm.time()

print(f"> Finished solving ODE for {T} simulated seconds, {P} data points in {e_time - s_time :.3f} real time seconds")

theta1 = sol[:, 0]
theta2 = sol[:, 1]

# Static Plots
plt.figure()
plt.plot(t, theta1, color='lime', label='θ₁')
plt.plot(t, theta2, color='red', label='θ₂')
plt.xlabel('Time (s)')
plt.ylabel('Angle (radians)')
plt.title('System Angles')
plt.legend()
plt.grid(True)
plt.show()

plt.figure()
plt.plot(theta1, sol[:, 2], color='lime', label='θ₁')
plt.plot(theta2, sol[:, 3], color='red', label='θ₂')
plt.title('Angular Velocities Plotted Against Their Angles')
plt.xlabel('Angle (radians)')
plt.ylabel('Angular Velocities')
plt.legend()
plt.show()

# Precompute Cartesian coordinates
X1 = l1 * np.sin(theta1)
Y1 = -l1 * np.cos(theta1)
X2 = X1 + l2 * np.sin(theta2)
Y2 = Y1 - l2 * np.cos(theta2)

# Static trace plot
plt.figure()
plt.scatter(0, 0, color='grey')
plt.plot(X1, Y1, color='cyan', linewidth=2, label='Mass 1 Path')
plt.plot(X2, Y2, color='purple', linewidth=2, label='Mass 2 Path')
plt.title("Paths followed by masses over time")
plt.xlabel('X Space')
plt.ylabel("Y Space")
plt.legend()
plt.grid(True)
plt.show()

# --- Animation Rendering using FuncAnimation ---
print("> Generating animation 'dpendu.gif'...")
s_time_2 = tm.time()

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-(l1 + l2 + 1), (l1 + l2 + 1))
ax.set_ylim(-(l1 + l2 + 1), 5)
ax.set_aspect('equal')
ax.grid(True)
ax.set_title("Double Pendulum Physical Simulation")

# Interactive plot elements
pivot = ax.scatter(0, 0, s=20, color='grey', zorder=5)
mass1 = ax.scatter([], [], s=80, color='green', zorder=5)
mass2 = ax.scatter([], [], s=80, color='red', zorder=5)
rod1, = ax.plot([], [], color='blue', linewidth=1.5)
rod2, = ax.plot([], [], color='orange', linewidth=1.5)
trail, = ax.plot([], [], color='purple', alpha=0.7, linewidth=1)

def init():
    rod1.set_data([], [])
    rod2.set_data([], [])
    mass1.set_offsets(np.empty((0, 2)))
    mass2.set_offsets(np.empty((0, 2)))
    trail.set_data([], [])
    return rod1, rod2, mass1, mass2, trail

def update(frame):
    x1, y1 = X1[frame], Y1[frame]
    x2, y2 = X2[frame], Y2[frame]
    
    rod1.set_data([0, x1], [0, y1])
    rod2.set_data([x1, x2], [y1, y2])
    
    mass1.set_offsets([[x1, y1]])
    mass2.set_offsets([[x2, y2]])
    
    trail.set_data(X2[:frame+1], Y2[:frame+1])
    
    return rod1, rod2, mass1, mass2, trail

anim = animation.FuncAnimation(
    fig, update, frames=len(t), init_func=init, interval=1000/fps, blit=True
)

output_file = "dpendu.gif"
anim.save(output_file, writer='pillow', fps=fps)
plt.close(fig)

e_time_2 = tm.time()
print(f"> Saved animation to '{output_file}' in {e_time_2 - s_time_2 :.3f} seconds.")
print("> End")