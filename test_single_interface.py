#!/usr/bin/env python3
"""
Test energy conservation for single interface (Air|Glass)
Check: P_toward + |P_away| = 2/3 in same layer
"""
import numpy as np
from te_integrate_plot import compute_power_flux_kp

# 2-layer: Air | Glass (single interface)
n_list = [1.0, 1.5]
d_list = [0.0]  # interface at z=0
wl = 0.65
k0 = 2 * np.pi / wl
ka = k0

print('Single interface (Air|Glass), source in air')
print('=' * 60)

for z_src in [-0.05, -0.1, -0.5, -1.0, -2.0]:
    # z_obs1: between source and interface (z > z_src)
    z_obs1 = z_src + 0.02
    # z_obs2: away from interface (z < z_src)
    z_obs2 = z_src - 0.02
    
    P1 = compute_power_flux_kp(n_list, d_list, 0, z_src, 0, z_obs1, k0)  # toward interface
    P2 = compute_power_flux_kp(n_list, d_list, 0, z_src, 0, z_obs2, k0)  # away from interface
    
    # P1 is positive (toward +z), P2 is negative (toward -z)
    # Total = P1 - P2 (since P2 points in -z direction)
    total = (P1 - P2) / ka
    
    print(f'z_src = {z_src:.2f} um:')
    print(f'  P_toward/ka = {P1/ka:.4f}, P_away/ka = {P2/ka:.4f}')
    print(f'  Total = P_toward - P_away = {total:.4f} (should be 2/3 = {2/3:.4f})')
    print()

print('\nNote: Total should be 2/3 if energy is conserved in source layer.')
print('      Some energy goes to glass, but our integral only covers kp < ka_air.')
