#!/usr/bin/env python3
"""
Demo 2: Supermassive Binary Black Hole -- test particle orbits
Restricted 3-body problem: ~10 test particles in the time-dependent
potential of two SMBHs on a Keplerian orbit.

Output: bh_xy.pdf, bh_xz.pdf, bh_yz.pdf
All projections use identical axis limits for fair comparison.
"""
import agama, numpy as np, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

agama.setUnits(length=1, velocity=1, mass=1)

def main():
    print("=== Demo 2: Binary BH test particle orbits ===")

    M_total = 1e8; q = 0.3; a = 0.005; ecc = 0.3
    T_orb = 2 * np.pi * np.sqrt(a**3 / (agama.G * M_total))

    pot_binary = agama.Potential(
        type='KeplerBinary', mass=M_total, binary_q=q,
        binary_sma=a, binary_ecc=ecc)
    pot_host = agama.Potential(type='Plummer', mass=1e10, scaleRadius=0.5)
    pot = agama.Potential(pot_host, pot_binary)

    np.random.seed(42)
    ics = []
    R_vals = [0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20]
    for R in R_vals:
        vc = np.sqrt(abs(pot.force([[R, 0, 0]])[0, 0] * R)) if R > 0 else 0
        phi = np.random.uniform(0, 2*np.pi)
        x = R * np.cos(phi); y = R * np.sin(phi); z = R * np.random.uniform(-0.1, 0.1)
        vx = -vc * np.sin(phi); vy = vc * np.cos(phi); vz = vc * np.random.uniform(-0.05, 0.05)
        ics.append([x, y, z, vx, vy, vz])
    ics = np.array(ics)
    print(f"  {len(ics)} test particles")

    T_int = 50 * T_orb
    print("  Integrating ...")
    times, traj = agama.orbit(potential=pot, ic=ics, time=T_int,
                               trajsize=5001, separateTime=True)

    # Determine common axis limits from all trajectories
    all_xyz = np.vstack([t for t in traj])
    x_min, x_max = np.min(all_xyz[:,0]), np.max(all_xyz[:,0])
    y_min, y_max = np.min(all_xyz[:,1]), np.max(all_xyz[:,1])
    z_min, z_max = np.min(all_xyz[:,2]), np.max(all_xyz[:,2])
    lim = max(abs(x_min), abs(x_max), abs(y_min), abs(y_max), abs(z_min), abs(z_max))

    colours = plt.cm.viridis(np.linspace(0, 1, len(ics)))

    for proj, fname in [([0, 1], 'bh_xy.pdf'), ([0, 2], 'bh_xz.pdf'), ([1, 2], 'bh_yz.pdf')]:
        fig, ax = plt.subplots(figsize=(7, 7))
        for i in range(len(ics)):
            ax.plot(traj[i][:, proj[0]], traj[i][:, proj[1]], '-',
                    lw=0.4, color=colours[i], alpha=0.7)
        ax.plot(0, 0, 'k+', ms=12, mew=2)
        ax.set_xlabel(['x [kpc]', 'x [kpc]', 'y [kpc]'][proj[0]])
        ax.set_ylabel(['y [kpc]', 'z [kpc]', 'z [kpc]'][proj[1]])
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_aspect('equal')
        ax.set_title(f'Binary BH: {fname[:-4]}')
        fig.tight_layout(); fig.savefig(fname); plt.close()
        print(f"  Saved {fname}")

    print("=== Demo 2 complete ===")

if __name__ == '__main__':
    main()
