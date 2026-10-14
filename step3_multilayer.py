"""
Step 3: Transient conduction through a multilayer wall with convective boundary conditions on both sides 

Uses a finite-volume formulation so that thermal conductivity can jump discontinuously at the material interface
 """


import numpy as np 
import matplotlib.pyplot as plt

#--- Layer definitions ------------
# (thickness [m], k [W/m-k], rho [kg.m^3], cp [J,kg-K])
#k = thermal conductivity, rho = density, cp = specific heat capacity
layer_A = dict(thickness = 0.01, k = 15, rho = 8000.0, cp = 500) #Specific heat capacity of Stainless Steel
layer_B = dict(thickness = .03, k =0.04, rho = 100, cp = 1400) #Specific heat capacity of wood/other insulators

# ---- Boundary Conditions
T_hot = 180 #deg C, hot process fluid
h_hot = 500 #W/m^2-K, hot-side convection coefficient
T_cold = 25 #deg C, ambient air
h_cold = 10 #cold-side convection coefficient
T_initial = 25.0 #wall starts at ambient temperature

#----- Numerical Definitions ---------------------------
nx_A = 20
nx_B = 60
L_A = layer_A["thickness"]
L_B = layer_B["thickness"]
dx = L_A / nx_A
assert abs(L_B / nx_B - dx)<1e-12 

nx = nx_A + nx_B
x = np.linspace (0, L_A + L_B, nx)

k_arr = np.where (x<L_A, layer_A["k"], layer_B ["k"])
rho_arr = np.where(x < L_A, layer_A["rho"], layer_B["rho"])
cp_arr = np.where(x < L_A, layer_A["cp"], layer_B["cp"])

# --- Timer Step --------
alpha_A = layer_A["k"] / (layer_A["rho"] * layer_A["cp"])
alpha_B = layer_B["k"] / (layer_B["rho"] * layer_B["cp"])
alpha_max = max(alpha_A,alpha_B)
dt_max = 0.5 * dx**2 / alpha_max
dt = .25*dt_max
t_final = 1800 #seconds
nt = int(t_final/dt)
print(f"alpha_A={alpha_A:.3e}, alpha_B={alpha_B:.3e} m^2/s")
print(f"dx={dx:.3e} m, dt={dt:.4f} s, steps={nt}")

def k_face (kL, kR):
    """Harmonic mean conductivity at the face between two control volumes."""
    return 2 * kL*kR/(kL+kR)

# --- Initialize --------
T = np.full(nx, T_initial)
snapshot_times = (0,60,300,900,1800)
snapshots = {}
snap_idx = 0

# --- Time funciton ---------
for n in range (nt + 1): 
    t = n*dt
    if snap_idx < len(snapshot_times) and t>= snapshot_times[snap_idx]:
        snapshots[snapshot_times[snap_idx]] = T.copy()
        snap_idx += 1

    T_new = T.copy()

    #Interior nodes
    for i in range (1,nx-1):
        k_e = k_face(k_arr[i], k_arr[i+1])
        k_w = k_face(k_arr[i-1],k_arr[i])
        q_net = k_e * (T[i+1] - T[i])/ dx - k_w * (T[i] - T[i-1]) / dx
        T_new[i] = T[i] + dt / (rho_arr[i] * cp_arr[i] * dx) * q_net

    #Left boundary 
    k_e = k_face(k_arr[0],k_arr[1])
    q_net = h_hot * (T_hot - T[0]) - k_e *(T[0] - T[1])/dx
    T_new[0] = T[0] + dt / (rho_arr[0] * cp_arr[0] * dx / 2) * q_net

    # Right boundary: convective to ambient, half control volume
    k_w = k_face(k_arr[-2], k_arr[-1])
    q_net = k_w * (T[-2] - T[-1]) / dx - h_cold * (T[-1] - T_cold)
    T_new[-1] = T[-1] + dt / (rho_arr[-1] * cp_arr[-1] * dx / 2) * q_net

    T = T_new

# ---- Plot ----------------------------------------------------------------------
plt.figure(figsize=(7.5, 4.8))
for t_snap, T_snap in snapshots.items():
    plt.plot(x * 1000, T_snap, label=f"t = {t_snap:>4d} s")
plt.axvline(L_A * 1000, color='gray', linestyle='--', linewidth=1, label="material interface")
plt.xlabel("Position (mm)")
plt.ylabel("Temperature (°C)")
plt.title("Transient Conduction Through a Composite Insulated Wall")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("step3_multilayer.png", dpi=150)
print("Saved plot to step3_multilayer.png")
print("Final surface temps: hot-side =", T[0], " cold-side =", T[-1])