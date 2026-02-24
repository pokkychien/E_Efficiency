#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

def r_te(qi, qj):
    """TE Fresnel reflection (using q = sqrt(kp**2 - eps*k0**2)/k0)."""
    return (qi - qj) / (qi + qj)

def r_tm(qi, qj, ei, ej):
    """TM Fresnel reflection (using q = sqrt(kp**2 - eps*k0**2)/k0)."""
    return (qi/ei - qj/ej) / (qi/ei + qj/ej)

def RF_multilayer(n_list, d_list, polarization, k0, kp, target_layer):
    e_list = [n**2 for n in n_list]
    kp2 = kp**2
    r_eff = 0.0 + 0.0j
    for L in range(len(n_list)-1, target_layer, -1):
        eps_up   = e_list[L-1]
        eps_down = e_list[L]
        q_up   = np.sqrt(kp2 - eps_up   * k0**2 + 0j) / k0
        q_down = np.sqrt(kp2 - eps_down * k0**2 + 0j) / k0
        xL     = np.exp(-q_up * k0 * d_list[L-1])
        if polarization.upper() == "TE":
            wn = q_down / q_up
        else:
            wn = (q_down / q_up) * (eps_up / eps_down)
        wm = 1.0
        wp = wm + wn
        wq = wm - wn
        al = (wp + r_eff * wq) / 2.0
        bl = (wq + r_eff * wp) / 2.0
        r_eff = xL * (bl / al) * xL
    return r_eff

def RB_multilayer(n_list, d_list, polarization, k0, kp, target_layer):
    e_list = [n**2 for n in n_list]
    kp2 = kp**2
    r_eff = 0.0 + 0.0j
    for L in range(0, target_layer):
        eps_up   = e_list[L]
        eps_down = e_list[L+1]
        q_up   = np.sqrt(kp2 - eps_up   * k0**2 + 0j) / k0
        q_down = np.sqrt(kp2 - eps_down * k0**2 + 0j) / k0
        xL     = np.exp(-q_up * k0 * d_list[L])
        if polarization.upper() == "TE":
            wn = q_up / q_down
        else:
            wn = (q_up / q_down) * (eps_down / eps_up)
        wm = 1.0
        wp = wm + wn
        wq = wm - wn
        bl = (wq + wp * xL * r_eff * xL)
        al = (wp + wq * xL * r_eff * xL)
        r_eff = bl / al
    return r_eff

def main():
    n_list = [1.0, 1.5, 1.0]
    d_list = [0.0, 500e-9]  # thicknesses (最後一層無厚度)
    wavelength = 533e-9  # 固定波長
    k0 = 2 * np.pi / wavelength
    # 掃 kp (從 0 到 n1*k0)
    points = 500
    kp = np.linspace(0, n_list[0]*k0, points)
    R_RF = np.empty(points, dtype=float)
    R_RB = np.empty(points, dtype=float)
    for i in range(points):
        R_RF[i] = np.abs(RF_multilayer(n_list, d_list, "TE", k0, kp[i], target_layer=1))**2
        R_RB[i] = np.abs(RB_multilayer(n_list, d_list, "TM", k0, kp[i], target_layer=1))**2
    plt.figure(figsize=(9,5))
    plt.plot(kp/k0, R_RF, label="RF from top (TE)")
    plt.plot(kp/k0, R_RB, label="RB from bottom (TM)")
    plt.xlabel(r"$k_\parallel/k_0$")
    plt.ylabel("Reflectance")
    plt.title(f"multi-layer reflectance (wavelength = {wavelength*1e9:.0f} nm)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
