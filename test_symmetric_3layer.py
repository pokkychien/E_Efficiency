#!/usr/bin/env python3
"""
三層對稱結構能量守恆與 RF/RB 對稱性測試
- n = [1.0, 1.5, 1.0]
- d = [0.0, 0.5]
- z_src = 0.25 (中間層正中央)
- z_obs = 0.24, 0.26 (上下各一點)
"""
import numpy as np
from te_greens import gyy_TE, dgyy_dzobs_TE
from tm_greens import gxx_TM, gzx_TM, dgxx_dzobs_TM

n_list = [1.0, 1.5, 1.0]
d_list = [0.0, 0.5]
z_src = 0.25
wl = 0.6
k0 = 2 * np.pi / wl

# 積分範圍
kp_min = 0.01 * k0
kp_max = k0 * n_list[1] * 0.98
num_k = 200
kp_vals = np.linspace(kp_min, kp_max, num_k)

from te_integrate_plot import compute_power_flux_kp

def flux(layer_obs, z_obs):
    # 這裡 layer_src 固定為 1, z_src 為全域變數
    return compute_power_flux_kp(
        n_list, d_list, 1, z_src, layer_obs, z_obs, k0
    )

print("=== Symmetric 3-layer structure ===")
print(f"n = {n_list}, d = {d_list}")
print(f"z_src = {z_src}")

# 檢查 source 上下方能量流
f_up = flux(1, z_src + 0.02)
f_down = flux(1, z_src - 0.02)
print(f"Flux above source (z={z_src+0.02}): {f_up:.6f}")
print(f"Flux below source (z={z_src-0.02}): {f_down:.6f}")
print(f"Sum (should be ~0): {f_up + f_down:.6e}")

# 檢查 RF/RB 絕對值對稱性
kp_test = k0 * 0.5
Gyy_up = gyy_TE(n_list, d_list, 1, z_src, 2, 0.5, k0, kp_test)
Gyy_down = gyy_TE(n_list, d_list, 1, z_src, 0, 0.0, k0, kp_test)
print(f"\nRF (up, abs): {abs(Gyy_up):.6f}")
print(f"RB (down, abs): {abs(Gyy_down):.6f}")
print(f"Difference: {abs(abs(Gyy_up) - abs(Gyy_down)):.2e}")
