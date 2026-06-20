#!/usr/bin/env python3
"""
Demo 4: Triaxial elliptical galaxy -- box orbits and loop (tube) orbits
All orbits are fully 3D with z-position and z-velocity.
Box orbits fill a volume with alternating sign of all three angular momentum components.
Loop (tube) orbits preserve sign of at least one angular momentum component.

Output: ellip_box_orbit.pdf, ellip_loop_orbit.pdf
All projections use identical axis limits.
"""
import agama, numpy as np, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

agama.setUnits(length=1, velocity=1, mass=1)
TU = 0.978

def main():
    print("=== Demo 4: Triaxial elliptical -- box & loop orbits ===")
    pot = agama.Potential(
        type='Dehnen', mass=1e11, scaleRadius=2.0, gamma=1.0,
        axisRatioY=0.8, axisRatioZ=0.6)
    print("  Triaxial Dehnen (1:0.8:0.6), M=1e11, a=2 kpc")

    # Generate many fully 3D ICs with z-velocity
    np.random.seed(456)
    candidates = []
    for R in [1.5, 2.0, 3.0, 4.0]:
        vc = (-R * pot.force([[R, 0, 0]])[0, 0])**0.5
        for f_tan in [0.15, 0.30, 0.50, 0.70, 0.90]:
            v_phi = vc * f_tan
            v_r = vc * 0.3 * (np.random.rand() - 0.5)
            v_z = vc * 0.3 * (np.random.rand() - 0.5)   # z-velocity!
            for phi in [0, np.pi/3, 2*np.pi/3]:
                for tilt in [0.1, 0.3, 0.5]:              # z-position
                    x = R * np.cos(phi)
                    y = R * np.sin(phi) * 0.8
                    z = R * tilt * 0.3
                    vx = v_r*np.cos(phi) - v_phi*np.sin(phi)
                    vy = v_r*np.sin(phi) + v_phi*np.cos(phi)
                    candidates.append([x, y, z, vx, vy, v_z])
    candidates = np.array(candidates)
    print(f"  Testing {len(candidates)} fully 3D ICs ...")

    T_int = 10.0 / TU
    _, trajs = agama.orbit(potential=pot, ic=candidates, time=T_int,
                           trajsize=3001, separateTime=True)

    # Classification: box vs loop based on angular momentum flips
    box_scores = []; loop_scores = []
    for i in range(len(candidates)):
        tr = trajs[i]
        L = np.cross(tr[:, :3], tr[:, 3:6])
        L_mag = np.mean(np.sqrt(np.sum(L**2, axis=1)))
        flips = sum(np.sum(np.diff(np.sign(L[:, k])) != 0) for k in range(3))
        # Box: many flips in all 3 components
        box_score = flips
        # Loop: few flips, non-zero angular momentum
        loop_score = -flips + L_mag * 0.01
        box_scores.append(box_score)
        loop_scores.append(loop_score)

    box_idx = np.argmax(box_scores)
    loop_idx = np.argmax(loop_scores)
    print(f"  Box candidate {box_idx}: box_score={box_scores[box_idx]:.0f}")
    print(f"  Loop candidate {loop_idx}: loop_score={loop_scores[loop_idx]:.0f}")

    # Common axis limits
    lim = max(np.max(np.abs(trajs[box_idx][:,:3])), np.max(np.abs(trajs[loop_idx][:,:3])))

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
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
            ax.set_aspect('equal'); ax.set_xlabel(xl); ax.set_ylabel(yl)
        axes[0].set_title(f'Triaxial: {label} Orbit (3D)')
        fig.tight_layout(); fig.savefig(fname); plt.close()
        print(f"  Saved {fname}")

    print("=== Demo 4 complete ===")

if __name__ == '__main__':
    main()
