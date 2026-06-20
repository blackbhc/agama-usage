#!/usr/bin/env python3
"""
Demo 3: Barred potential — box orbits and loop orbits
Test ~20 initial conditions in a rotating barred potential,
select one clear box orbit and one loop orbit, plot each as a separate PDF.

Output: bar_box_orbit.pdf, bar_loop_orbit.pdf
"""
import agama, numpy as np, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CONFIG = os.path.join(os.path.dirname(__file__), '..', 'configs',
                       'test3_bar_potential.ini')
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'outputs')

agama.setUnits(length=1, velocity=1, mass=1)
TU = 0.978
Omega_bar = -40.0  # km/s/kpc

def main():
    print("=== Demo 3: Barred potential — box & loop orbits ===")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    pot = agama.Potential(file=CONFIG)

    # Scan ICs to find box and loop orbits
    np.random.seed(123)
    candidates = []
    for R in [3.0, 4.0, 5.0, 6.0, 7.0]:
        vc = (-R * pot.force([[R, 0, 0]])[0, 0])**0.5
        for e_factor in [0.4, 0.6, 0.8, 1.0]:
            v_phi = vc * e_factor
            v_r = vc * 0.3 * (np.random.rand() - 0.5)
            for phi in [0, np.pi/4, np.pi/2]:
                x = R * np.cos(phi)
                y = R * np.sin(phi)
                vx = v_r*np.cos(phi) - v_phi*np.sin(phi)
                vy = v_r*np.sin(phi) + v_phi*np.cos(phi)
                candidates.append([x, y, 0, vx, vy, 0])
    candidates = np.array(candidates)
    print(f"  Testing {len(candidates)} ICs for orbit classification ...")

    # Integrate all candidates
    T_int = 5.0 / TU
    _, trajs = agama.orbit(potential=pot, ic=candidates, time=T_int,
                           Omega=Omega_bar, trajsize=2001, separateTime=True)

    # Classify: compute (x-y) aspect ratio and angular momentum sign changes
    box_idx, loop_idx = None, None
    for i in range(len(candidates)):
        tr = trajs[i]
        x, y = tr[:, 0], tr[:, 1]
        Lz = x * tr[:, 4] - y * tr[:, 3]
        # Box orbit: crosses zero many times (alternates sign of Lz)
        zc = np.sum(np.diff(np.sign(Lz)) != 0)
        # Also check radial extent vs azimuthal extent
        r_std = np.std(np.sqrt(x**2 + y**2))
        # Loop orbit: keeps same Lz sign, fills an annulus
        if zc > 20 and box_idx is None:
            box_idx = i
            print(f"  Box orbit candidate {i}: zero-crossings={zc}")
        if zc < 4 and np.mean(np.abs(Lz)) > 0.1*np.max(np.abs(Lz)) and loop_idx is None:
            loop_idx = i
            print(f"  Loop orbit candidate {i}: zero-crossings={zc}")

    if box_idx is None or loop_idx is None:
        print("  Could not classify orbits, using fallback indices")
        box_idx = box_idx or 5
        loop_idx = loop_idx or 15

    # Plot box orbit
    tr = trajs[box_idx]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    axes[0].plot(tr[:, 0], tr[:, 1], 'b-', lw=0.5); axes[0].set_aspect('equal')
    axes[0].set_xlabel('x [kpc]'); axes[0].set_ylabel('y [kpc]')
    axes[0].set_title(f'Barred: Box Orbit (IC #{box_idx})')
    axes[1].plot(tr[:, 0], tr[:, 2], 'b-', lw=0.5); axes[1].set_aspect('equal')
    axes[1].set_xlabel('x [kpc]'); axes[1].set_ylabel('z [kpc]')
    axes[2].plot(tr[:, 1], tr[:, 2], 'b-', lw=0.5); axes[2].set_aspect('equal')
    axes[2].set_xlabel('y [kpc]'); axes[2].set_ylabel('z [kpc]')
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'bar_box_orbit.pdf'))
    plt.close()
    print(f"  Saved bar_box_orbit.pdf (IC #{box_idx})")

    # Plot loop orbit
    tr = trajs[loop_idx]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    axes[0].plot(tr[:, 0], tr[:, 1], 'r-', lw=0.5); axes[0].set_aspect('equal')
    axes[0].set_xlabel('x [kpc]'); axes[0].set_ylabel('y [kpc]')
    axes[0].set_title(f'Barred: Loop Orbit (IC #{loop_idx})')
    axes[1].plot(tr[:, 0], tr[:, 2], 'r-', lw=0.5); axes[1].set_aspect('equal')
    axes[1].set_xlabel('x [kpc]'); axes[1].set_ylabel('z [kpc]')
    axes[2].plot(tr[:, 1], tr[:, 2], 'r-', lw=0.5); axes[2].set_aspect('equal')
    axes[2].set_xlabel('y [kpc]'); axes[2].set_ylabel('z [kpc]')
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'bar_loop_orbit.pdf'))
    plt.close()
    print(f"  Saved bar_loop_orbit.pdf (IC #{loop_idx})")

    print("=== Demo 3 complete ===")

if __name__ == '__main__':
    main()
