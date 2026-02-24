
from __future__ import annotations

import numpy as np
from scipy import special
from scipy.integrate import quad
import matplotlib.pyplot as plt

from te_greens import gyy_TE, dgyy_dzobs_TE
from tm_greens import gxx_TM, dgxx_dzobs_TM, gzx_TM


def compute_power_flux_kp(
    n_list, d_list,
    layer_src: int, z_src: float,
    layer_obs: int, z_obs: float,
    k0: float, limit: int = 200,
) -> float:
    """
    Compute power flux by integrating over k_parallel from 0 to ka.
    Uses scipy.integrate.quad for adaptive integration.
    
    P = Im{ ∫_0^{ka} k_parallel * [
        i * k_parallel * G_xx * conj(G_zx) 
        + G_xx * conj(dG_xx/dz) 
        + G_yy * conj(dG_yy/dz)
    ] dk_parallel }
    
    where ka = k0 * n_obs (observation layer)
    
    Returns the imaginary part (real power flux).
    """
    n_obs = n_list[layer_obs]
    ka = k0 * n_obs  # integration upper limit
    
    def integrand(kp):
        gyy = gyy_TE(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
        dgyy = dgyy_dzobs_TE(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
        gxx = gxx_TM(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
        dgxx = dgxx_dzobs_TM(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
        gzx = gzx_TM(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
        val = kp * (
            1j * kp * gxx * np.conj(gzx)
            + gxx * np.conj(dgxx)
            + gyy * np.conj(dgyy)
        )
        return np.imag(val)
    P, _ = quad(integrand, 0, 0.995*ka, limit=limit)
    
    return P


def test_free_space_scaling():
    """
    Test power flux in free space (RF = RB = 0) for different ka values.
    Check if P is proportional to ka.
    """
    # Free space: single layer with n=1
    # Use uniform medium to ensure RF = RB = 0
    
    wl = 0.650  # wavelength in um
    k0 = 2 * np.pi / wl
    
    layer_src, z_src = 0, 0.0
    layer_obs, z_obs = 0, 0.1  # small separation to avoid cusp
    
    # Test different refractive indices (different ka = k0 * n)
    n_values = [1.0, 1.5, 2.0, 2.5, 3.0]
    results = []
    
    print("Free space power flux scaling test:")
    print("=" * 50)
    print(f"k0 = {k0:.4f}, z_src = {z_src}, z_obs = {z_obs}")
    print("-" * 50)
    print(f"{'n':>6} {'ka':>10} {'P':>15} {'P/ka':>15}")
    print("-" * 50)
    
    for n in n_values:
        n_list = [n, n, n, n]  # uniform medium
        d_list = [0.0, 2.0, 1.5]
        
        ka = k0 * n
        P = compute_power_flux_kp(
            n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0
        )
        
        results.append((n, ka, P))
        print(f"{n:>6.2f} {ka:>10.4f} {P:>15.6e} {P/ka:>15.6e}")
    
    print("-" * 50)
    print("\nIf P ∝ ka, then P/ka should be constant (should be 1/3).")
    
    return results

