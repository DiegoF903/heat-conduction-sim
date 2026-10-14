""" 
Step 2: Validate the explicit STCS solver against an exact analytical solution for a rod initially at uniform temperature T_i using Fourier series.
"""

import numpy as np
import matplotlib.pyplot as plt

# --- Physical parameters --------
L= 0.05
alpha = 1.1e-4
T_s = 25.0
T_i = 100.0

#--- Numerical Definitions ---------------
nx = 41
dx = L / (nx-1)
x = np.linspace(0,L,nx)

dt_max = 0.5 * dx**2 / alpha 
dt = 0.4* dt_max 
t_check = 3.0
nt = int(t_check / dt)

r = alpha * dt / dx**2

# ---- FTCS Solution ------------
T = np.full(nx, T_i)
T[0] = T_s
T[-1] = T_s
for n in range(nt):
    T_new = T.copy()
    T_new[1:-1] = T[1:-1] + r * (T[2:] - 2 * T[1:-1] + T[:-2])
    T_new[0] = T_s
    T_new[-1] = T_s
    T = T_new

t_actual = nt * dt

#---Fourier Solution, first 200 terms -----------
def analytical(x,t,L,alpha,T_i,T_s,n_terms=200):
    T = np.full_like(x,T_s,dtype = float)
    for k in range (1,n_terms,+1):
        n= 2*k-1  #odd terms of the sin
        coeff = 4.0 / (n*np.pi)
        T += (T_i - T_s) * coeff * np.sin (n*np.pi*x/L) * np.exp (-(n*np.pi/L) **2*alpha*t)
    return T

T_exact = analytical(x,t_actual,L, alpha, T_i,T_s)

# ---Comparing analytical and fourier solutions ----------------
error = np.abs(T - T_exact)
print(f"Comparing at t = {t_actual:.4f} s")
print(f"Max absolute error: {error.max():.4f} deg C")
print(f"Mean absolute error: {error.mean():.4f} deg C")

plt.figure(figsize=(7, 4.5))
plt.plot(x * 1000, T_exact, 'k-', linewidth=2, label="Analytical (Fourier series)")
plt.plot(x * 1000, T, 'o', markersize=4, color='tab:red', label="Numerical (FTCS)")
plt.xlabel("Position (mm)")
plt.ylabel("Temperature (°C)")
plt.title(f"Numerical vs. Analytical Solution at t = {t_actual:.2f} s")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("step2_validation.png", dpi=150)
print("Saved plot to step2_validation.png")
