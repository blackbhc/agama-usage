#!/usr/bin/env python3
"""
Demo 4: Elliptical galaxy -- box orbits and loop (tube) orbits
Triaxial Dehnen potential representing an elliptical galaxy.
Test various ICs, select one box orbit and one loop orbit.

Output: ellip_box_orbit.pdf, ellip_loop_orbit.pdf (saved in current directory)
"""
import agama, numpy as np, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

agama.setUnits(length=1, velocity=1, mass=1)
TU = 0.978

def main():
    print("=== Demo 4: Elliptical galaxy -- box & loop orbits ===")
    pot = agama.Potential(
        type='Dehnen', mass=1e11, scaleRadius=2.0, gamma=1.0,
        axisRatioY=0.8, axisRatioZ=0.6)
    print("  Triaxial Dehnen (1:0.8:0.6), M=1e11, a=2 kpc")

    np.random.seed(456)
    candidates = []
    for R in [1.5, 2.0, 3.0, 4.0]:
        vc = (-R * pot.force([[R, 0, 0]])[0, 0])**0.5
        for e_factor in [0.3, 0.5, 0.7, 0.9]:
            v_phi = vc * e_factor
            for phi in [0, np.pi/3, 2*np.pi/3]:
                for tilt in [0, np.pi/6]:
                    x = R * np.cos(phi); y = R * np.sin(phi) * 0.8
                    z = R * np.sin(tilt) * 0.3
                    vx = -v_phi * np.sin(phi); vy = v_phi * np.cos(phi) * 0.8
                    vz = v_phi * np.sin(tilt) * 0.2
                    candidates.append([x, y, z, vx, vy, vz])
    candidates = np.array(candidates)
    print(f"  Testing {len(candidates)} ICs ...")

    T_int = 8.0 / TU
    _, trajs = agama.orbit(potential=pot, ic=candidates, time=T_int,
                           trajsize=2001, separateTime=True)

    box_idx = loop_idx = None
    for i in range(len(candidates)):
        tr = trajs[i]
        L = np.cross(tr[:, :3], tr[:, 3:6])
        L_mag = np.mean(np.sqrt(np.sum(L**2, axis=1)))
        flips = sum(np.sum(np.diff(np.sign(L[:, k])) != 0) for k in range(3))
        if flips > 60 and box_idx is None:
            box_idx = i; print(f"  Box candidate {i}: {flips} flips")
        if flips < 20 and L_mag > 0.1 and loop_idx is None:
            loop_idx = i; print(f"  Loop candidate {i}: {flips} flips")

    if box_idx is None: box_idx = 0
    if loop_idx is None: loop_idx = 2

    for idx, label, colour, fname in [
        (box_idx, 'Box', 'blue', 'ellip_box_orbit.pdf'),
        (loop_idx, 'Loop/Tube', 'red', 'ellip_loop_orbit.pdf')]:
        tr = trajs[idx]
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
        for ax, (i, j), xl, yl in [
            (axes[0], (0, 1), 'x [kpc]', 'y [kpc]'),
            (axes[1], (0, 2), 'x [kpc]', 'z [kpc]'),
            (axes[2], (1, 2), 'y [kpc]', 'z [kpc]')]:
            ax.plot(tr[:, i], tr[:, j], '-', lw=0.5, color=colour)
            ax.set_aspect('equal'); ax.set_xlabel(xl); ax.set_ylabel(yl)
        axes[0].set_title(f'Triaxial: {label} Orbit')
        fig.tight_layout(); fig.savefig(fname); plt.close()
        print(f"  Saved {fname}")

    print("=== Demo 4 complete ===")

if __name__ == '__main__':
    main()
