#!/usr/bin/env python3
"""
Demo 3: Barred potential -- box orbits and loop orbits
Test ~20 initial conditions in a rotating barred potential,
select one clear box orbit and one loop orbit.

Output: bar_box_orbit.pdf, bar_loop_orbit.pdf (saved in current directory)
"""
import agama, numpy as np, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CONFIG = os.path.join(os.path.dirname(__file__), 'bar_potential.ini')
agama.setUnits(length=1, velocity=1, mass=1)
TU = 0.978
Omega_bar = -40.0

def main():
    print("=== Demo 3: Barred potential -- box & loop orbits ===")
    pot = agama.Potential(file=CONFIG)
    np.random.seed(123)
    candidates = []
    for R in [3.0, 4.0, 5.0, 6.0, 7.0]:
        vc = (-R * pot.force([[R, 0, 0]])[0, 0])**0.5
        for e_factor in [0.4, 0.6, 0.8, 1.0]:
            v_phi = vc * e_factor
            v_r = vc * 0.3 * (np.random.rand() - 0.5)
            for phi in [0, np.pi/4, np.pi/2]:
                x = R * np.cos(phi); y = R * np.sin(phi)
                vx = v_r*np.cos(phi) - v_phi*np.sin(phi)
                vy = v_r*np.sin(phi) + v_phi*np.cos(phi)
                candidates.append([x, y, 0, vx, vy, 0])
    candidates = np.array(candidates)
    print(f"  Testing {len(candidates)} ICs ...")

    T_int = 5.0 / TU
    _, trajs = agama.orbit(potential=pot, ic=candidates, time=T_int,
                           Omega=Omega_bar, trajsize=2001, separateTime=True)

    box_idx = loop_idx = None
    for i in range(len(candidates)):
        tr = trajs[i]
        Lz = tr[:, 0] * tr[:, 4] - tr[:, 1] * tr[:, 3]
        zc = np.sum(np.diff(np.sign(Lz)) != 0)
        if zc > 20 and box_idx is None:
            box_idx = i; print(f"  Box candidate {i}: {zc} zero-crossings")
        if zc < 4 and np.mean(np.abs(Lz)) > 1 and loop_idx is None:
            loop_idx = i; print(f"  Loop candidate {i}: {zc} zero-crossings")

    if box_idx is None: box_idx = 5
    if loop_idx is None: loop_idx = 0

    for idx, label, colour, fname in [
        (box_idx, 'Box', 'blue', 'bar_box_orbit.pdf'),
        (loop_idx, 'Loop', 'red', 'bar_loop_orbit.pdf')]:
        tr = trajs[idx]
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
        for ax, (i, j), xl, yl in [
            (axes[0], (0, 1), 'x [kpc]', 'y [kpc]'),
            (axes[1], (0, 2), 'x [kpc]', 'z [kpc]'),
            (axes[2], (1, 2), 'y [kpc]', 'z [kpc]')]:
            ax.plot(tr[:, i], tr[:, j], '-', lw=0.5, color=colour)
            ax.set_aspect('equal'); ax.set_xlabel(xl); ax.set_ylabel(yl)
        axes[0].set_title(f'Barred: {label} Orbit')
        fig.tight_layout(); fig.savefig(fname); plt.close()
        print(f"  Saved {fname}")

    print("=== Demo 3 complete ===")

if __name__ == '__main__':
    main()
