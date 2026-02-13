"""
Bessel transform: G(k_parallel) -> G(rho)

Integrals:
  I1  = (1/π) ∫ k * J0(k*ρ) * G_yy^TE(k) dk
  I2  = (1/π) ∫ k * J0(k*ρ) * conj(dG_yy^TE/dz_obs(k)) dk
  I3  = (1/π) ∫ k * J2(k*ρ) * G_yy^TE(k) dk
  I4  = (1/π) ∫ k * J2(k*ρ) * conj(dG_yy^TE/dz_obs(k)) dk
  I5  = (1/π) ∫ k * J0(k*ρ) * G_xx^TM(k) dk
  I6  = (1/π) ∫ k * J0(k*ρ) * conj(dG_xx^TM/dz_obs(k)) dk
  I7  = (1/π) ∫ k * J2(k*ρ) * G_xx^TM(k) dk
  I8  = (1/π) ∫ k * J2(k*ρ) * conj(dG_xx^TM/dz_obs(k)) dk
  I9  = (1/π) ∫ k² * J0(k*ρ) * conj(G_zx^TM(k)) dk
  I10 = (1/π) ∫ k² * J2(k*ρ) * conj(G_zx^TM(k)) dk

Using scipy Bessel functions and split integration to avoid k=k0 singularity.
"""

from __future__ import annotations

import numpy as np
from scipy import special
import matplotlib.pyplot as plt

from te_greens import gyy_TE, dgyy_dzobs_TE
from tm_greens import gxx_TM, dgxx_dzobs_TM, gzx_TM


def compute_all_integrals(
    n_list, d_list,
    layer_src: int, z_src: float,
    layer_obs: int, z_obs: float,
    k0: float, rho: float,
    k_parallel_max: float, num_k: int,
    delta: float = 0.01,
) -> dict:
    """
    Compute all 10 integrals (I1-I10) using split integration.
    
    Split integration to avoid singularity at k_parallel = k0:
      - Region 1: [0, k0*(1-delta)]           — num_k points
      - Region 2: [k0*(1+delta), k_parallel_max] — num_k points
    
    delta: fractional gap around k0 (e.g. 0.01 means skip k0±1%)
    Each region gets num_k points independently.
    
    Returns dict with keys: I1, I2, ..., I10
    """
    
    k_left  = k0 * (1.0 - delta)   # end of region 1
    k_right = k0 * (1.0 + delta)   # start of region 2
    
    # Region 1: [0, k_left] with num_k points
    dk1 = k_left / num_k
    kps1 = (np.arange(num_k) + 0.5) * dk1
    
    # Region 2: [k_right, k_parallel_max] with num_k points
    dk2 = (k_parallel_max - k_right) / num_k
    kps2 = k_right + (np.arange(num_k) + 0.5) * dk2
    
    # Combine
    kps = np.concatenate([kps1, kps2])
    n_total = 2 * num_k
    
    # ========== Compute Green's functions ==========
    # G_yy^TE
    Gyykp = np.array([
        gyy_TE(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
        for kp in kps
    ], dtype=np.complex128)
    
    # dG_yy^TE/dz_obs
    DGyykp = np.array([
        dgyy_dzobs_TE(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
        for kp in kps
    ], dtype=np.complex128)
    
    # G_xx^TM
    Gxxkp = np.array([
        gxx_TM(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
        for kp in kps
    ], dtype=np.complex128)
    
    # dG_xx^TM/dz_obs  
    DGxxkp = np.array([
        dgxx_dzobs_TM(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
        for kp in kps
    ], dtype=np.complex128)
    
    # G_zx^TM
    Gzxkp = np.array([
        gzx_TM(n_list, d_list, layer_src, z_src, layer_obs, z_obs, k0, kp)
        for kp in kps
    ], dtype=np.complex128)
    
    # ========== Bessel functions (scipy) ==========
    J0 = special.jv(0, kps * rho)
    J2 = special.jv(2, kps * rho)
    
    # ========== Build integrands ==========
    integrand_1  = kps * J0 * Gyykp
    integrand_2  = kps * J0 * np.conj(DGyykp)
    integrand_3  = kps * J2 * Gyykp
    integrand_4  = kps * J2 * np.conj(DGyykp)
    integrand_5  = kps * J0 * Gxxkp
    integrand_6  = kps * J0 * np.conj(DGxxkp)
    integrand_7  = kps * J2 * Gxxkp
    integrand_8  = kps * J2 * np.conj(DGxxkp)
    integrand_9  = kps**2 * J0 * np.conj(Gzxkp)
    integrand_10 = kps**2 * J2 * np.conj(Gzxkp)
    
    # ========== Integrate (split regions) ==========
    pref = 1.0 / np.pi
    
    def split_trapz(integrand):
        """Integrate over both regions and sum."""
        I_r1 = np.trapz(integrand[:num_k], kps[:num_k])
        I_r2 = np.trapz(integrand[num_k:], kps[num_k:])
        return pref * (I_r1 + I_r2)
    
    return {
        "I1":  split_trapz(integrand_1),
        "I2":  split_trapz(integrand_2),
        "I3":  split_trapz(integrand_3),
        "I4":  split_trapz(integrand_4),
        "I5":  split_trapz(integrand_5),
        "I6":  split_trapz(integrand_6),
        "I7":  split_trapz(integrand_7),
        "I8":  split_trapz(integrand_8),
        "I9":  split_trapz(integrand_9),
        "I10": split_trapz(integrand_10),
    }


def gyy_TE_rho(
    n_list, d_list,
    layer_src: int, z_src: float,
    layer_obs: int, z_obs: float,
    k0: float, rho: float,
    k_parallel_max: float, num_k: int,
    delta: float = 0.01,
) -> complex:
    """
    Compute I1 only (Bessel-Fourier transform of G_yy^TE).
    For backward compatibility.
    """
    integrals = compute_all_integrals(
        n_list, d_list, layer_src, z_src, layer_obs, z_obs,
        k0, rho, k_parallel_max, num_k, delta
    )
    return integrals["I1"]


def full_efficiency_integrand(
    n_list, d_list,
    layer_src: int, z_src: float,
    layer_obs: int, z_obs: float,
    k0: float, rho: float,
    k_parallel_max: float, num_k: int,
    delta: float = 0.01,
) -> complex:
    """
    Compute the full efficiency integrand combining all terms:
      I1*I2 + I3*I4 + I5*I6 + I7*I8 
      + 1j*I5*I9 + 1j*I5*I10 
      + I5*I2 + I7*I4 + I1*I6 + I3*I8 
      + 1j*I1*I9 + 1j*I1*I10
    """
    I = compute_all_integrals(
        n_list, d_list, layer_src, z_src, layer_obs, z_obs,
        k0, rho, k_parallel_max, num_k, delta
    )
    
    return (
        I["I1"] * I["I2"] + I["I3"] * I["I4"] +
        I["I5"] * I["I6"] + I["I7"] * I["I8"] +
        1j * I["I5"] * I["I9"] + 1j * I["I5"] * I["I10"] +
        I["I5"] * I["I2"] + I["I7"] * I["I4"] +
        I["I1"] * I["I6"] + I["I3"] * I["I8"] +
        1j * I["I1"] * I["I9"] + 1j * I["I1"] * I["I10"]
    )


def demo():
    """
    Demo plot showing I1 (G_yy^TE Bessel transform).
    """
    # Geometry
    n_list = [1.0, 1.0, 1.0, 1.0]
    d_list = [0.0, 2.0, 1.5]
    
    layer_src, z_src = 0, -0.35
    layer_obs, z_obs = 0, -0.20
    
    # k0
    wl = 0.650
    k0 = 2 * np.pi / wl
    
    # integration
    k_parallel_max = 10.0 * k0
    num_k = 1001
    delta = 0.0001
    
    # rho sweep
    rhos = np.linspace(0.0, 20.0, 500)
    G_rho = np.array([
        gyy_TE_rho(n_list, d_list, layer_src, z_src, layer_obs, z_obs,
                   k0, r, k_parallel_max, num_k, delta=delta)
        for r in rhos
    ])

    # plot
    plt.figure()
    plt.plot(rhos, np.real(G_rho), label="Re I1")
    plt.plot(rhos, np.imag(G_rho), label="Im I1")
    plt.axhline(0, color='k', ls='--', lw=0.8)
    plt.xlabel(r"$\rho$ (µm)")
    plt.ylabel(r"$I_1(\rho)$")
    plt.legend()
    plt.title("I1: Bessel transform of $G_{yy}^{TE}$")
    plt.tight_layout()
    plt.show()


def demo_all_integrals():
    """
    Demo plot showing all integrals I1-I10.
    """
    # Geometry
    n_list = [1.0, 1.0, 1.0, 1.0]
    d_list = [0.0, 2.0, 1.5]
    
    layer_src, z_src = 0, -0.35
    layer_obs, z_obs = 0, -0.20
    
    # k0
    wl = 0.650
    k0 = 2 * np.pi / wl
    
    # integration
    k_parallel_max = 1.0 * k0
    num_k = 501
    delta = 0.0001
    
    # rho sweep
    rhos = np.linspace(0.01, 10.0, 200)
    
    all_results = [
        compute_all_integrals(
            n_list, d_list, layer_src, z_src, layer_obs, z_obs,
            k0, r, k_parallel_max, num_k, delta=delta
        )
        for r in rhos
    ]
    
    # Extract each integral
    labels = ["I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9", "I10"]
    
    fig, axes = plt.subplots(2, 5, figsize=(15, 6))
    axes = axes.flatten()
    
    for i, label in enumerate(labels):
        data = np.array([r[label] for r in all_results])
        axes[i].plot(rhos, np.real(data), label="Re")
        axes[i].plot(rhos, np.imag(data), label="Im")
        axes[i].axhline(0, color='k', ls='--', lw=0.5)
        axes[i].set_title(label)
        axes[i].set_xlabel(r"$\rho$ (µm)")
        axes[i].legend(fontsize=8)
    
    plt.tight_layout()
    plt.show()


def demo_full_efficiency():
    """
    Demo plot showing the full efficiency integrand.
    """
    # Geometry
    n_list = [1.0, 1.0, 1.0, 1.0]
    d_list = [0.0, 2.0, 1.5]
    
    layer_src, z_src = 0, -0.30
    layer_obs, z_obs = 0, -0.20
    
    # k0
    wl = 0.650
    k0 = 2 * np.pi / wl
    
    # integration
    k_parallel_max = 10.0 * k0
    num_k = 501
    delta = 0.0001
    
    # rho sweep
    rhos = np.linspace(0.01, 5.0, 400)
    
    eff = np.array([
        full_efficiency_integrand(
            n_list, d_list, layer_src, z_src, layer_obs, z_obs,
            k0, r, k_parallel_max, num_k, delta=delta
        )
        for r in rhos
    ])
    
    # Cylindrical coordinate integration: 2π ∫ eff(ρ) * ρ dρ
    integrand_cyl = eff * rhos
    total_integral = 2 * np.pi * np.trapz(integrand_cyl, rhos)
    
    plt.figure()
    plt.plot(rhos, np.real(eff), label="Re")
    plt.plot(rhos, np.imag(eff), label="Im")
    plt.axhline(0, color='k', ls='--', lw=0.8)
    plt.xlabel(r"$\rho$ (µm)")
    plt.ylabel("Full Efficiency Integrand")
    plt.legend()
    plt.title("Full Efficiency Integrand (all terms)")
    
    # Display the integrated value on the plot
    text_str = (f"$2\\pi \\int_0^{{\\rho_{{max}}}} f(\\rho) \\rho \\, d\\rho$\n"
                f"= {np.real(total_integral):.6f} + {np.imag(total_integral):.6f}i")
    plt.text(0.95, 0.95, text_str, transform=plt.gca().transAxes,
             fontsize=10, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    print(f"Cylindrical integral: {total_integral}")
    print(f"  Real part: {np.real(total_integral)}")
    print(f"  Imag part: {np.imag(total_integral)}")
    
    plt.tight_layout()
    plt.show()


def plot_Grho_as_2D(G_rho, rhos, extent_um=3.5, N=401, which="real"):
    """
    Make a 2D map from radial data G(ρ) by coordinate transform + interpolation.

    extent_um: plot x,y in [-extent_um, +extent_um] (µm)
    N: grid points per axis
    which: "real" or "imag" or "abs"
    """
    x = np.linspace(-extent_um, extent_um, N)
    y = np.linspace(-extent_um, extent_um, N)
    X, Y = np.meshgrid(x, y, indexing="xy")
    R = np.sqrt(X**2 + Y**2)

    if which == "real":
        g1d = np.real(G_rho)
        label = "Re G"
    elif which == "imag":
        g1d = np.imag(G_rho)
        label = "Im G"
    elif which == "abs":
        g1d = np.abs(G_rho)
        label = "|G|"
    else:
        raise ValueError('which must be "real", "imag", or "abs"')

    G2 = np.interp(R.ravel(), rhos, g1d, left=0.0, right=0.0).reshape(R.shape)

    plt.figure()
    im = plt.imshow(
        G2,
        origin="lower",
        extent=[-extent_um, extent_um, -extent_um, extent_um],
        aspect="equal",
    )
    plt.colorbar(im, label=label)
    plt.xlabel("x (µm)")
    plt.ylabel("y (µm)")
    plt.title(f"{label} from radial G(ρ)")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Choose which demo to run:
    # demo()            # I1 only (original behavior)
    # demo_all_integrals()   # All 10 integrals
     demo_full_efficiency() # Full efficiency integrand
