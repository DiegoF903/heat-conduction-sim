""" 
Step1: Explicit finite-difference solver for 1D transient heat conduction in a uniform rod with fixed temperature ends.

Governing equation: dT/dt = alpha * d2T/dx2
"""

import os
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import matplotlib.pyplot as plt

#---Physical Definitions  --------------
L = 0.05            #Fixed rod length, m
alpha = 1.1e-4      #thermal diffusivity of copper, m^2/s
T_left = 100.0      # held left end temperature, deg C
T_right = 25.0      # held constant right end temperature, deg C
T_initial = 25.0    # initial rod temperature, deg C

#---Numerical definitions -------------------------------------------------

nx = 41 # number of spatial nodes
dx = L / (nx - 1)
x = np.linspace (0, L, nx) 

# Stability limit for explicit FTCS (forward time difference used to approximate the first time derivative and an approximate second spatial derivative) in 1D: alpha*dt/dx^2 <= 0.5
dt_max = 0.5 * dx**2 / alpha
dt = 0.4 * dt_max     # 80% of the limit
t_final = 8.0         # seconds
nt = int(t_final / dt)

r = alpha * dt / dx**2
print(f"dx = {dx: 0.3e} m, dt = {dt: 0.3e} s, r = alpha*dt/dx^2 = {r: 0.3f} (must be <= 0.5)")
print (f"Number of time steps: {nt}")

# --- Initial Condition --------------------
T = np.full(nx, T_initial)
T[0] = T_left
T[-1] = T_right

# ----- Explicit FTCS ----------
frames = [T.copy()]
snapshots = {}
snapshots_times = [0, 0.5, 1, 2, 4, 8]
snap_idx = 0

for n in range(nt + 1):
    t = n * dt
    if snap_idx < len(snapshots_times) and t >= snapshots_times[snap_idx]:
        snapshots[snapshots_times[snap_idx]] = T.copy()
        snap_idx += 1

    T_new = T.copy()
    T_new[1:-1] = T[1:-1] + r * (T[2:] - 2 * T[1:-1] + T[:-2])
    T_new[0] = T_left
    T_new[-1] = T_right
    T = T_new
    frames.append(T.copy())
    
    

# --- Plot ----------------
plt.figure(figsize=(7, 4.5))
for t_snap, T_snap in snapshots.items():
    plt.plot(x * 1000, T_snap, marker='o', markersize=3, label=f"t = {t_snap:.1f} s")
plt.xlabel("Position (mm)")
plt.ylabel("Temperature (°C)")
plt.title("1D Transient Conduction in a Rod (Explicit FTCS)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("step1_explicit.png", dpi=150)
print("Saved plot to step1_explicit.png")
print("Final temperature profile (first/mid/last node):", T[0], T[nx//2], T[-1])

# --- GIF Animation ----------------
fig, ax = plt.subplots(figsize=(7, 4.5))


print("x shape:", x.shape)
print("frames length:", len(frames))
print("frames[0] shape:", np.shape(frames[0]))
line, = ax.plot(x * 1000, frames[0], 'o-', markersize=3)

ax.set_xlabel("Position (mm)")
ax.set_ylabel("Temperature (°C)")
ax.set_title("1D Transient Conduction in a Rod (Explicit FTCS)")
ax.set_xlim(0, L * 1000)
ax.set_ylim(T_initial - 5, T_left + 5)
ax.grid(alpha=0.3)

time_text = ax.text(
    0.02, 0.95, "", transform=ax.transAxes, va="top"
)

def update(frame):
    line.set_ydata(frames[frame])
    time_text.set_text(f"Time = {frame * dt:.3f} s")
    return line, time_text

ani = FuncAnimation(
    fig,
    update,
    frames=len(frames),
    interval=50,
    blit=False,
    repeat=True
)

output_file = os.path.abspath("step1_animation.gif")

ani.save(
    output_file,
    writer=PillowWriter(fps=20)
)

plt.close(fig)

print(f"GIF saved to: {output_file}")
print(f"File size: {os.path.getsize(output_file) / 1024:.1f} KB")