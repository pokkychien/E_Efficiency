#!/usr/bin/env python3
"""
Scan z_obs position for fixed source, plot |P/ka| vs z_obs
Three-layer structure for complex RF
"""
import numpy as np
import matplotlib.pyplot as plt
from te_integrate_plot import compute_power_flux_kp

# Parameters
wl = 0.65  # um (650 nm)
k0 = 2 * np.pi / wl

# Three-layer structure: Air | Glass film (200nm) | Air
n_list = [1.0, 1.5, 1.0]  # air | glass | air
d_list = [0.0, 0.2]  # glass film thickness = 200 nm

# Fixed source in air (layer 0)
layer_src = 0
z_src = -0.1  # 100 nm above interface

# Observation layer (same as source, in air)
layer_obs = 0
ka = k0 * n_list[layer_obs]

# Scan z_obs from -0.12 to -2.0 um
z_obs_arr = np.linspace(-0.12, -2.0, 50)

print("Three-layer structure: Air | Glass (200nm) | Air")
print("Scanning z_obs (Source in AIR, fixed at z_src = -0.1 um)...")
print("=" * 60)

P_ka_arr = []
for z_obs in z_obs_arr:
    P = compute_power_flux_kp(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0)
    P_ka = np.abs(P / ka)
    P_ka_arr.append(P_ka)
    print(f"  z_obs = {z_obs:.3f} um, |P/ka| = {P_ka:.4f}")

P_ka_arr = np.array(P_ka_arr)

# Plot
plt.figure(figsize=(8, 5))
plt.plot(z_obs_arr, P_ka_arr, 'b-', linewidth=2)
plt.axhline(y=1/3, color='r', linestyle='--', label='Free space (1/3)')
plt.xlabel('z_obs (um)', fontsize=12)
plt.ylabel('|P/ka|', fontsize=12)
plt.title(f'Power flux vs observation position\n(λ={wl*1000:.0f}nm, z_src={z_src}um, Air|Glass(200nm)|Air)', fontsize=12)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('power_vs_zobs.png', dpi=150)
print("\nPlot saved as power_vs_zobs.png")
plt.show()
