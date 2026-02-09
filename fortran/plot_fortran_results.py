"""
Plot Fortran results and optionally compare with Python results.
"""
import numpy as np
import matplotlib.pyplot as plt

def load_fortran_data(filename='gyy_rho_fortran.dat'):
    """Load data from Fortran output file."""
    data = np.loadtxt(filename, comments='#')
    rho = data[:, 0]
    re_g = data[:, 1]
    im_g = data[:, 2]
    return rho, re_g + 1j * im_g

def plot_fortran_only(filename='gyy_rho_fortran.dat'):
    """Plot Fortran results."""
    rho, G = load_fortran_data(filename)
    
    plt.figure(figsize=(10, 6))
    plt.plot(rho, np.real(G), 'b-', label='Re(G) Fortran', linewidth=1.5)
    plt.plot(rho, np.imag(G), 'r-', label='Im(G) Fortran', linewidth=1.5)
    plt.axhline(0, color='k', linestyle='--', linewidth=0.5)
    plt.xlabel(r'$\rho$ (µm)')
    plt.ylabel(r'$G(\rho)$')
    plt.title('Fortran: Bessel-integrated Green\'s function')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('gyy_rho_fortran.png', dpi=150)
    plt.show()

def compare_with_python():
    """Compare Fortran results with Python calculation."""
    import sys
    sys.path.insert(0, '..')
    
    from te_integrate_plot import gyy_TE_rho
    
    # Load Fortran data
    rho_f, G_f = load_fortran_data('gyy_rho_fortran.dat')
    
    # Python parameters (same as Fortran)
    n_list = [1.0, 1.0, 1.0, 1.0]
    d_list = [0.0, 2.0, 1.5]
    layer_src = 0  # Python 0-based
    layer_obs = 0
    z_src = -0.35
    z_obs = -0.20
    wl = 0.650
    k0 = 2.0 * np.pi / wl
    k_parallel_max = 3.5 * k0
    num_k = 2001
    
    # Compute Python results
    print("Computing Python results for comparison...")
    G_p = np.array([
        gyy_TE_rho(
            n_list, d_list,
            layer_src, z_src,
            layer_obs, z_obs,
            k0,
            rho=float(r),
            k_parallel_max=k_parallel_max,
            num_k=num_k,
            use_trapz=True,
        )
        for r in rho_f
    ], dtype=np.complex128)
    
    # Plot comparison
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Real part
    axes[0, 0].plot(rho_f, np.real(G_f), 'b-', label='Fortran', linewidth=1.5)
    axes[0, 0].plot(rho_f, np.real(G_p), 'r--', label='Python', linewidth=1.5)
    axes[0, 0].set_xlabel(r'$\rho$ (µm)')
    axes[0, 0].set_ylabel(r'Re(G)')
    axes[0, 0].set_title('Real Part')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Imaginary part
    axes[0, 1].plot(rho_f, np.imag(G_f), 'b-', label='Fortran', linewidth=1.5)
    axes[0, 1].plot(rho_f, np.imag(G_p), 'r--', label='Python', linewidth=1.5)
    axes[0, 1].set_xlabel(r'$\rho$ (µm)')
    axes[0, 1].set_ylabel(r'Im(G)')
    axes[0, 1].set_title('Imaginary Part')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Difference (Real)
    diff_re = np.real(G_f) - np.real(G_p)
    axes[1, 0].plot(rho_f, diff_re, 'g-', linewidth=1)
    axes[1, 0].set_xlabel(r'$\rho$ (µm)')
    axes[1, 0].set_ylabel(r'$\Delta$ Re(G)')
    axes[1, 0].set_title(f'Real Part Difference (max = {np.max(np.abs(diff_re)):.2e})')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Difference (Imag)
    diff_im = np.imag(G_f) - np.imag(G_p)
    axes[1, 1].plot(rho_f, diff_im, 'g-', linewidth=1)
    axes[1, 1].set_xlabel(r'$\rho$ (µm)')
    axes[1, 1].set_ylabel(r'$\Delta$ Im(G)')
    axes[1, 1].set_title(f'Imag Part Difference (max = {np.max(np.abs(diff_im)):.2e})')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.suptitle('Fortran vs Python Comparison')
    plt.tight_layout()
    plt.savefig('fortran_vs_python.png', dpi=150)
    plt.show()
    
    print(f"\nMax |Re diff|: {np.max(np.abs(diff_re)):.6e}")
    print(f"Max |Im diff|: {np.max(np.abs(diff_im)):.6e}")
    print(f"Max |G_f|: {np.max(np.abs(G_f)):.6e}")
    print(f"Relative error: {np.max(np.abs(G_f - G_p)) / np.max(np.abs(G_f)):.6e}")

if __name__ == "__main__":
    import os
    
    # Change to fortran directory
    fortran_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(fortran_dir)
    
    print("Plotting Fortran results...")
    plot_fortran_only()
    
    # Uncomment to compare with Python:
    # compare_with_python()
