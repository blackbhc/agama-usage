#!/usr/bin/env python3
"""
Demo 3: Barred potential -- box orbits and loop orbits in the corotating frame.
The bar is elongated along the x-axis. In the frame rotating with the bar,
box orbits oscillate along the bar's major axis, while loop orbits circulate.

Output: bar_box_orbit.pdf, bar_loop_orbit.pdf
All projections use identical axis limits.
"""
import agama, numpy as np, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CONFIG = os.path.join(os.path.dirname(__file__), 'bar_potential.ini')
agama.setUnits(length=1, velocity=1, mass=1)
TU = 0.978
Omega_bar = -40.0  # km/s/kpc, bar pattern speed (corotating frame)

def main():
    print("=== Demo 3: Barred potential -- box & loop orbits ===")
    pot = agama.Potential(file=CONFIG)

    # Generate a range of ICs: vary radius, tangential fraction, azimuth
    # Box orbits (corotating frame): radial motion dominant, little rotation
    # Loop orbits: significant tangential motion, net Lz
    candidates = []
    info = []
    for R in [3.0, 4.0, 5.0, 6.0]:
        vc = (-R * pot.force([[R, 0, 0]])[0, 0])**0.5
        for f_tan in [0.05, 0.15, 0.30, 0.50, 0.70, 0.90]:  # tangential fraction of vc
            v_phi = vc * f_tan
            v_r = vc * 0.4 * (np.random.rand() - 0.5)
            for phi in [0, np.pi/4, np.pi/3, np.pi/2]:
                x = R * np.cos(phi); y = R * np.sin(phi)
                vx = v_r*np.cos(phi) - v_phi*np.sin(phi)
                vy = v_r*np.sin(phi) + v_phi*np.cos(phi)
                candidates.append([x, y, 0, vx, vy, 0])
                info.append((R, f_tan, phi))
    candidates = np.array(candidates)
    print(f"  Testing {len(candidates)} ICs ...")

    T_int = 6.0 / TU
    _, trajs = agama.orbit(potential=pot, ic=candidates, time=T_int,
                           Omega=Omega_bar, trajsize=3001, separateTime=True)

    # Classification in the corotating frame:
    # Box orbit: Lz zero-crossings many times, x-extent ~ y-extent or x > y
    # Loop orbit: Lz keeps sign, shape is more circular
    box_scores = []   # higher = better box
    loop_scores = []  # higher = better loop
    for i in range(len(candidates)):
        tr = trajs[i]
        x, y = tr[:, 0], tr[:, 1]
        Lz = x * tr[:, 4] - y * tr[:, 3]
        zc = np.sum(np.diff(np.sign(Lz)) != 0)
        # Mean Lz magnitude
        mean_Lz = np.mean(np.abs(Lz))
        # Axial ratio in x-y: box orbits along bar (x) should have x_std > y_std
        x_std, y_std = np.std(x), np.std(y)
        axial_ratio = min(x_std, y_std) / max(x_std, y_std) if max(x_std, y_std) > 0 else 1
        # Box score: many zero-crossings, axial ratio not too small (filled)
        box_score = zc * (1 - axial_ratio) if zc > 10 else 0
        # Loop score: few zero-crossings, large mean Lz, decent axial ratio (roundish)
        loop_score = -zc + mean_Lz / 100 if zc < 10 else 0
        box_scores.append(box_score)
        loop_scores.append(loop_score)

    box_idx = np.argmax(box_scores)
    loop_idx = np.argmax(loop_scores)
    print(f"  Box orbit:  IC #{box_idx}  (R={info[box_idx][0]:.1f}, "
          f"f_tan={info[box_idx][1]:.2f}, zc={np.sum(np.diff(np.sign(trajs[box_idx][:,0]*trajs[box_idx][:,4]-trajs[box_idx][:,1]*trajs[box_idx][:,3]))!=0)})")
    print(f"  Loop orbit: IC #{loop_idx} (R={info[loop_idx][0]:.1f}, "
          f"f_tan={info[loop_idx][1]:.2f}, zc={np.sum(np.diff(np.sign(trajs[loop_idx][:,0]*trajs[loop_idx][:,4]-trajs[loop_idx][:,1]*trajs[loop_idx][:,3]))!=0)})")

    # Determine common axis limits
    box_tr = trajs[box_idx]; loop_tr = trajs[loop_idx]
    lim = max(np.max(np.abs(box_tr[:,:3])), np.max(np.abs(loop_tr[:,:3])))

    for idx, label, colour, fname in [
        (box_idx, 'Box (along bar axis)', 'blue', 'bar_box_orbit.pdf'),
        (loop_idx, 'Loop', 'red', 'bar_loop_orbit.pdf')]:
        tr = trajs[idx]
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
        for ax, (i, j), xl, yl in [
            (axes[0], (0, 1), 'x [kpc]', 'y [kpc]'),
            (axes[1], (0, 2), 'x [kpc]', 'z [kpc]'),
            (axes[2], (1, 2), 'y [kpc]', 'z [kpc]')]:
            ax.plot(tr[:, i], tr[:, j], '-', lw=0.5, color=colour)
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
            ax.set_aspect('equal'); ax.set_xlabel(xl); ax.set_ylabel(yl)
        axes[0].set_title(f'Barred (corotating): {label}')
        fig.tight_layout(); fig.savefig(fname); plt.close()
        print(f"  Saved {fname}")

    print("=== Demo 3 complete ===")

if __name__ == '__main__':
    main()
