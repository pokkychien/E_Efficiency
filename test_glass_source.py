import numpy as np
from te_greens import gyy_TE, dgyy_dzobs_TE
from tm_greens import gxx_TM, dgxx_dzobs_TM, gzx_TM
from scipy.integrate import quad

# Two-layer structure:
# layer_0: air (n=1.0), z < 0
# layer_1: glass (n=1.5), z > 0

n_list = [1.0, 1.5]
d_list = [0.0]
wl = 0.650
k0 = 2 * np.pi / wl

print('=' * 70)
print('TEST 1: Source in AIR (n=1), away from interface')
print('        Reflection coefficient from air->glass is small (~0.04)')
print('=' * 70)

layer_src = 0  # air
layer_obs = 0  # air
z_src = -0.1   # in air, 100 nm from interface
z_obs = -0.12  # 20 nm away from interface

ka = k0 * n_list[layer_obs]

def integrand(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp):
    gyy = gyy_TE(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
    dgyy = dgyy_dzobs_TE(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
    gxx = gxx_TM(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
    dgxx = dgxx_dzobs_TM(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
    gzx = gzx_TM(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
    val = kp * (1j * kp * gxx * np.conj(gzx) + gxx * np.conj(dgxx) + gyy * np.conj(dgyy))
    return np.imag(val)

# With interface
P1, _ = quad(lambda kp: integrand(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp), 0, ka, limit=200)
# Free space
P1_free, _ = quad(lambda kp: integrand([1.0], [], 0, z_src, 0, z_obs, k0, kp), 0, ka, limit=200)

print(f'Source at z_src = {z_src}, Obs at z_obs = {z_obs}')
print(f'With interface: P/ka = {P1/ka:.6f}')
print(f'Free space:     P/ka = {P1_free/ka:.6f}')
print(f'|P/ka| with interface: {abs(P1/ka):.6f}')
print()

print('=' * 70)
print('TEST 2: Source in GLASS (n=1.5), away from interface')
print('        Reflection from glass->air is STRONGER')
print('        Plus total internal reflection for large angles!')
print('=' * 70)

layer_src = 1  # glass
layer_obs = 1  # glass
z_src = 0.1    # in glass, 100 nm from interface
z_obs = 0.12   # 20 nm further into glass (away from interface)

ka_glass = k0 * n_list[layer_obs]

# With interface
P2, _ = quad(lambda kp: integrand(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp), 0, ka_glass, limit=200)
# Free space (in glass)
P2_free, _ = quad(lambda kp: integrand([1.5], [], 0, z_src, 0, z_obs, k0, kp), 0, ka_glass, limit=200)

print(f'Source at z_src = {z_src}, Obs at z_obs = {z_obs}')
print(f'With interface: P/ka = {P2/ka_glass:.6f}')
print(f'Free space:     P/ka = {P2_free/ka_glass:.6f}')
print(f'|P/ka| with interface: {abs(P2/ka_glass):.6f}')
print()

print('=' * 70)
print('TEST 3: Source close to interface, observation far away')
print('        Maximizes reflected wave contribution')
print('=' * 70)

# Source in glass, very close to interface
layer_src = 1
layer_obs = 1
z_src = 0.02   # 20 nm from interface
z_obs = 0.2    # 200 nm into glass

P3, _ = quad(lambda kp: integrand(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp), 0, ka_glass, limit=200)
P3_free, _ = quad(lambda kp: integrand([1.5], [], 0, z_src, 0, z_obs, k0, kp), 0, ka_glass, limit=200)

print(f'Source at z_src = {z_src}, Obs at z_obs = {z_obs}')
print(f'With interface: P/ka = {P3/ka_glass:.6f}')
print(f'Free space:     P/ka = {P3_free/ka_glass:.6f}')
print(f'|P/ka| with interface: {abs(P3/ka_glass):.6f}')
print()

print('=' * 70)
print('SUMMARY')
print('=' * 70)
print('If |P/ka| > 1/3, reflection enhances power')
print('If |P/ka| < 1/3, interference reduces power')
print()
print('Normal incidence Fresnel reflection coefficients:')
r_TE = (1.0 - 1.5) / (1.0 + 1.5)
r_TM = (1.5 - 1.0) / (1.5 + 1.0)  # glass->air perspective
print(f'  Air->Glass: r = {r_TE:.4f}, R = {r_TE**2:.4f}')
print(f'  Glass->Air: r = {r_TM:.4f}, R = {r_TM**2:.4f}')
print(f'  Critical angle for TIR (glass->air): {np.arcsin(1.0/1.5)*180/np.pi:.1f} degrees')
