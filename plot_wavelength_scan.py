import numpy as np
import matplotlib.pyplot as plt
from te_integrate_plot import compute_power_flux_kp

# Three-layer structure: Air | Glass (500nm) | Air  (Fabry-Perot)
n_list = [1.0, 1.50, 1.0]
d_list = [0.0, 0.5]  # glass film thickness = 500 nm

# Source in AIR, observation away from interface
layer_src = 1
layer_obs = 2
z_src = 0.25   # 250 nm from interface (far field)
z_obs = 3.0  # 1500 nm behind source (away from interface)

# Wavelength scan (finer resolution for FP resonance)
wavelengths = np.linspace(0.3, 1.0, 500)  # 300 nm to 1000 nm

P_ka_values = []
free_space_values = []

print('Three-layer: Air | Glass (500nm) | Air')
print('Scanning wavelength...')
for wl in wavelengths:
    k0 = 2 * np.pi / wl
    ka = k0 * n_list[layer_obs]
    
    P = compute_power_flux_kp(n_list, d_list, layer_src, z_src, layer_obs, z_obs, ka)
    P_ka_values.append(abs(P / ka))
    free_space_values.append(1/3)
    print(f'  wl = {wl:.3f} um, |P/ka| = {abs(P/ka):.4f}')

# Plot
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(wavelengths * 1000, P_ka_values, 'b-', linewidth=1.5, label='Air|Glass(500nm)|Air')
ax.axhline(y=1/3, color='r', linestyle='--', linewidth=2, label='Free space (1/3)')
ax.set_xlabel('Wavelength (nm)', fontsize=12)
ax.set_ylabel('|P/ka|', fontsize=12)
ax.set_title(f'Fabry-Perot: Power flux vs Wavelength\nAir|Glass(500nm)|Air, z_src={z_src*1000:.0f}nm', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim([300, 1000])
ax.set_ylim([0, 0.6])

plt.tight_layout()
plt.savefig('power_vs_wavelength.png', dpi=150)
plt.show()

print(f'\nPlot saved as power_vs_wavelength.png')
