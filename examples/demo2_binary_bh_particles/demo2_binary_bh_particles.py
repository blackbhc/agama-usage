#!/usr/bin/env python3
"""
Demo 2: Supermassive Binary Black Hole -- test particle orbits
Restricted 3-body problem: ~10 test particles in the time-dependent
potential of two SMBHs on a Keplerian orbit.

Output: bh_xy.pdf, bh_xz.pdf, bh_yz.pdf (saved in current directory)
"""
import agama, numpy as np, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

agama.setUnits(length=1, velocity=1, mass=1)
TU = 0.978

def main():
    print("=== Demo 2: Binary BH test particle orbits ===")

    M_total = 1e8; q = 0.3; a = 0.005; ecc = 0.3
    T_orb = 2 * np.pi * np.sqrt(a**3 / (agama.G * M_total))
    print(f"  Binary: M={M_total:.1e}, q={q}, a={a:.4f} kpc, e={ecc}")

    pot_binary = agama.Potential(
        type='KeplerBinary', mass=M_total, binary_q=q,
        binary_sma=a, binary_ecc=ecc)
    pot_host = agama.Potential(type='Plummer', mass=1e10, scaleRadius=0.5)
    pot = agama.Potential(pot_host, pot_binary)

    np.random.seed(42)
    ics = []
    R_vals = [0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20]
    for R in R_vals:
        vc = np.sqrt(R * abs(pot.force([[R, 0, 0]])[0, 0]) / R) if R > 0 else 0
        phi = np.random.uniform(0, 2*np.pi)
        x = R * np.cos(phi); y = R * np.sin(phi); z = R * np.random.uniform(-0.1, 0.1)
        vx = -vc * np.sin(phi); vy = vc * np.cos(phi); vz = vc * np.random.uniform(-0.05, 0.05)
        ics.append([x, y, z, vx, vy, vz])
    ics = np.array(ics)
    print(f"  {len(ics)} test particles")

    T_int = 50 * T_orb
    print(f"  Integrating ...")
    times, traj = agama.orbit(potential=pot, ic=ics, time=T_int,
                               trajsize=5001, separateTime=True)

    colours = plt.cm.viridis(np.linspace(0, 1, len(ics)))

    for proj, fname in [([0, 1], 'bh_xy.pdf'), ([0, 2], 'bh_xz.pdf'), ([1, 2], 'bh_yz.pdf')]:
        fig, ax = plt.subplots(figsize=(7, 7))
        ax.set_aspect('equal')
        for i in range(len(ics)):
            ax.plot(traj[i][:, proj[0]], traj[i][:, proj[1]], '-', lw=0.4, color=colours[i], alpha=0.7)
        ax.plot(0, 0, 'k+', ms=12, mew=2)
        ax.set_xlabel(['x', 'x', 'y'][proj[0]] + ' [kpc]')
        ax.set_ylabel(['y', 'z', 'z'][proj[1]] + ' [kpc]')
        ax.set_title(f'Binary BH: {fname[:-4]}')
        fig.tight_layout()
        fig.savefig(fname)
        plt.close()
        print(f"  Saved {fname}")

    print("=== Demo 2 complete ===")

if __name__ == '__main__':
    main()
