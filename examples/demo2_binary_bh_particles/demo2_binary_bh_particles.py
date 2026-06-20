#!/usr/bin/env python3
"""
Demo 2: Supermassive Binary Black Hole — test particle orbits
Restricted 3-body problem: ~10 test particles in the time-dependent
potential of two SMBHs on a Keplerian orbit.

Output: bh_xy.pdf, bh_xz.pdf, bh_yz.pdf
"""
import agama, numpy as np, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'outputs')

agama.setUnits(length=1, velocity=1, mass=1)
TU = 0.978  # Gyr / time-unit

def main():
    print("=== Demo 2: Binary BH test particle orbits ===")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Binary SMBH parameters
    M_total = 1e8      # Msun
    q = 0.3            # mass ratio m2/m1
    a = 0.005          # semi-major axis (kpc) ~10 pc
    ecc = 0.3          # eccentricity
    T_orb = 2 * np.pi * np.sqrt(a**3 / (agama.G * M_total))
    print(f"  Binary: M={M_total:.1e}, q={q}, a={a:.4f} kpc, "
          f"e={ecc}, T={T_orb:.3f} time-units ({T_orb*TU:.3f} Gyr)")

    # Binary BH potential (time-dependent Kepler orbit)
    pot_binary = agama.Potential(
        type='KeplerBinary',
        mass=M_total, binary_q=q,
        binary_sma=a, binary_ecc=ecc)

    # Combine with a shallow Plummer core for the host galaxy
    pot_host = agama.Potential(type='Plummer', mass=1e10, scaleRadius=0.5)
    pot = agama.Potential(pot_host, pot_binary)

    # 10 test particles on near-circular orbits at various radii
    np.random.seed(42)
    ics = []
    R_vals = [0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20]
    for R in R_vals:
        vc = np.sqrt(R * pot.force(np.array([[R, 0, 0]]))[0, 0] / R) if R > 0 else 0
        phi = np.random.uniform(0, 2*np.pi)
        inc = np.random.uniform(-0.2, 0.2)
        x = R * np.cos(phi)
        y = R * np.sin(phi)
        z = R * inc * 0.3
        vx = -vc * np.sin(phi)
        vy =  vc * np.cos(phi)
        vz = vc * inc * 0.1
        ics.append([x, y, z, vx, vy, vz])
    ics = np.array(ics)
    print(f"  {len(ics)} test particles initialised")

    # Integrate for 50 binary orbits
    T_int = 50 * T_orb
    print(f"  Integrating for {T_int:.2f} time-units ({T_int*TU:.1f} Gyr) ...")
    times, traj = agama.orbit(potential=pot, ic=ics, time=T_int,
                               trajsize=5001, separateTime=True)
    print("  Done.")

    # Colour each particle by initial radius
    colours = plt.cm.viridis(np.linspace(0, 1, len(ics)))

    # Figure 1: x-y
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect('equal')
    for i in range(len(ics)):
        tr = traj[i]
        ax.plot(tr[:, 0], tr[:, 1], '-', lw=0.4, color=colours[i],
                alpha=0.7, label=f'R0={R_vals[i]:.3f}' if i in [0,9] else '')
    ax.plot(0, 0, 'k+', ms=12, mew=2)
    ax.set_xlabel('x [kpc]'); ax.set_ylabel('y [kpc]')
    ax.set_title('Binary BH: x-y projection')
    ax.legend(fontsize=7, loc='upper right')
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'bh_xy.pdf'))
    plt.close()
    print("  Saved bh_xy.pdf")

    # Figure 2: x-z
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect('equal')
    for i in range(len(ics)):
        tr = traj[i]
        ax.plot(tr[:, 0], tr[:, 2], '-', lw=0.4, color=colours[i], alpha=0.7)
    ax.plot(0, 0, 'k+', ms=12, mew=2)
    ax.set_xlabel('x [kpc]'); ax.set_ylabel('z [kpc]')
    ax.set_title('Binary BH: x-z projection')
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'bh_xz.pdf'))
    plt.close()
    print("  Saved bh_xz.pdf")

    # Figure 3: y-z
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect('equal')
    for i in range(len(ics)):
        tr = traj[i]
        ax.plot(tr[:, 1], tr[:, 2], '-', lw=0.4, color=colours[i], alpha=0.7)
    ax.plot(0, 0, 'k+', ms=12, mew=2)
    ax.set_xlabel('y [kpc]'); ax.set_ylabel('z [kpc]')
    ax.set_title('Binary BH: y-z projection')
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'bh_yz.pdf'))
    plt.close()
    print("  Saved bh_yz.pdf")

    print("=== Demo 2 complete ===")

if __name__ == '__main__':
    main()
