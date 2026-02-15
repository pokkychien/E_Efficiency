"""
Scan source position closer to interface to show evanescent wave coupling loss.

When source approaches interface, evanescent waves (kp > ka) couple energy
to glass, causing total power in air to drop below 2/3.
"""

import numpy as np
import matplotlib.pyplot as plt
from te_integrate_plot import compute_power_flux_kp

# Single interface: Air | Glass
n_list = [1.0, 1.5]
d_list = [0.0]  # interface at z=0

wl = 0.650  # um
k0 = 2 * np.pi / wl
ka = k0 * 1.0  # ka in air

# Scan z_src from far (-2 um) to near (-0.02 um)
z_src_list = np.linspace(-2.0, -0.02, 50)

P_toward_list = []
P_away_list = []
total_list = []

print("Scanning source position toward interface...")
print("=" * 60)

for z_src in z_src_list:
    # Power toward interface (z_obs slightly above z_src)
    z_obs_toward = z_src + 0.02
    if z_obs_toward > -0.001:
        z_obs_toward = -0.001  # stay in air
    
    # Power away from interface (z_obs below z_src)
    z_obs_away = z_src - 0.02
    
    P_toward = compute_power_flux_kp(n_list, d_list, 0, z_src, 0, z_obs_toward, k0)
    P_away = compute_power_flux_kp(n_list, d_list, 0, z_src, 0, z_obs_away, k0)
    
    P_toward_list.append(P_toward / ka)
    P_away_list.append(P_away / ka)
    total_list.append((P_toward - P_away) / ka)

P_toward_list = np.array(P_toward_list)
P_away_list = np.array(P_away_list)
total_list = np.array(total_list)

# Print some key values
print(f"z_src = {z_src_list[0]:.2f} um (far): Total = {total_list[0]:.4f}")
print(f"z_src = {z_src_list[-1]:.2f} um (near): Total = {total_list[-1]:.4f}")
print(f"Missing energy at near field: {2/3 - total_list[-1]:.4f}")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: Individual powers
ax1 = axes[0]
ax1.plot(-z_src_list, P_toward_list, 'b-', label=r'$P_{\mathrm{toward}}/k_a$')
ax1.plot(-z_src_list, -P_away_list, 'r-', label=r'$|P_{\mathrm{away}}|/k_a$')
ax1.axhline(1/3, color='gray', linestyle='--', alpha=0.7, label='Free space (1/3)')
ax1.set_xlabel('Distance from interface (um)')
ax1.set_ylabel('P / ka')
ax1.set_title('Power in each direction')
ax1.legend()
ax1.set_xlim([0, 2])
ax1.grid(True, alpha=0.3)

# Right: Total power
ax2 = axes[1]
ax2.plot(-z_src_list, total_list, 'k-', linewidth=2, label='Total power in air')
ax2.axhline(2/3, color='green', linestyle='--', linewidth=2, label='Free space limit (2/3)')
ax2.fill_between(-z_src_list, total_list, 2/3, alpha=0.3, color='red', 
                  label='Evanescent loss to glass')
ax2.set_xlabel('Distance from interface (um)')
ax2.set_ylabel(r'$(P_{\mathrm{toward}} - P_{\mathrm{away}}) / k_a$')
ax2.set_title('Total power in air vs source position')
ax2.legend()
ax2.set_xlim([0, 2])
ax2.set_ylim([0.4, 0.75])
ax2.grid(True, alpha=0.3)

plt.suptitle('Evanescent wave coupling: Source near Air|Glass interface', fontsize=14)
plt.tight_layout()
plt.savefig('evanescent_loss.png', dpi=150)
plt.show()

print("\nPlot saved to evanescent_loss.png")
