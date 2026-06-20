#!/usr/bin/env python3
"""
Demo 3: Barred potential — box and loop orbits in the corotating frame.
agama.orbit(Omega=Ω) integrates in the rotating frame; returned (x,y,z)
are corotating coordinates. The static Ferrers bar (axisRatioY=0.35) is
aligned with the x-axis and stationary in these coordinates.

True box orbits: Lz oscillates through zero frequently AND orbit is
elongated along the bar (x_std >> y_std).
Loop orbits: Lz preserves sign, circulate around center.

Output: bar_box_orbit.pdf, bar_loop_orbit.pdf
"""
import agama, numpy as np, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

CONFIG = os.path.join(os.path.dirname(__file__), 'bar_potential.ini')
agama.setUnits(length=1, velocity=1, mass=1)
TU = 0.978
Omega_bar = -40.0  # km/s/kpc

# Bar structural parameters (must match bar_potential.ini)
BAR_AXIS_X = 4.0       # scaleRadius (kpc)
BAR_AXIS_RATIO = 0.35  # axisRatioY — y/x

def main():
    print("=== Demo 3: Barred potential — box & loop orbits (corotating) ===")
    pot = agama.Potential(file=CONFIG)

    # ----- generate ICs -----
    candidates = []; info = []

    # Box candidates: launched from bar major axis, nearly radial
    for R in [2.5, 3.0, 3.5, 4.0, 4.5]:
        vc = (-R * pot.force([[R, 0, 0]])[0, 0])**0.5
        for f_tan in [0.008, 0.024, 0.04, 0.08, 0.24]:
            for f_rad in [0.01, 0.02, 0.05, 0.10, 0.20]:
                v_phi = vc * f_tan
                v_r = vc * f_rad * (1 if np.random.rand() > 0.5 else -1)
                candidates.append([R, 0, 0, v_r, v_phi, 0])
                info.append(('box', R, f_tan, f_rad))

    # Loop candidates: from various positions with significant tangential velocity
    for R in [3.0, 4.0, 5.0, 6.0]:
        vc = (-R * pot.force([[R, 0, 0]])[0, 0])**0.5
        for f_tan in [0.60, 0.70, 0.80, 0.90]:
            for phi in [0, np.pi/6, np.pi/4, np.pi/3]:
                v_phi = vc * f_tan
                v_r = vc * 0.2 * (np.random.rand() - 0.5)
                x = R * np.cos(phi); y = R * np.sin(phi)
                vx = v_r*np.cos(phi) - v_phi*np.sin(phi)
                vy = v_r*np.sin(phi) + v_phi*np.cos(phi)
                candidates.append([x, y, 0, vx, vy, 0])
                info.append(('loop', R, f_tan, 0))

    candidates = np.array(candidates)
    print(f"  {len(candidates)} ICs, Omega={Omega_bar} km/s/kpc")

    # ----- integrate in corotating frame -----
    T_int = 4.0 / TU
    _, trajs = agama.orbit(potential=pot, ic=candidates, time=T_int,
                           Omega=Omega_bar, trajsize=5001, separateTime=True)

    # ----- CLASSIFICATION (corrected) -----
    # True box orbit in a barred potential:
    #   1. Lz sign changes frequently (zc large) — crosses center
    #   2. Elongated along bar major axis (x_std >> y_std)
    # Loop orbit:
    #   1. Lz sign rarely or never changes (zc small)
    #   2. Substantial |Lz| (circulates)
    box_score = []; loop_score = []
    for i in range(len(candidates)):
        tr = trajs[i]
        x, y = tr[:, 0], tr[:, 1]
        Lz = x * tr[:, 4] - y * tr[:, 3]
        zc = np.sum(np.diff(np.sign(Lz)) != 0)
        x_std, y_std = np.std(x), np.std(y)
        elongation = x_std / y_std if y_std > 0 else 1.0  # >1 means stretched along x
        mean_Lz = np.mean(np.abs(Lz))

        # Box: frequently crosses center AND stretched along bar
        is_box = zc > 40 and elongation > 1.5
        bs = zc * elongation if is_box else 0

        # Loop: steady circulation
        is_loop = zc < 5 and mean_Lz > 5
        ls = mean_Lz / (zc + 1) if is_loop else 0

        box_score.append(bs); loop_score.append(ls)

    box_idx = np.argmax(box_score)
    loop_idx = np.argmax(loop_score)
    print(f"  Box:  IC #{box_idx}  R0={info[box_idx][1]:.1f}  "
          f"f_tan={info[box_idx][2]:.3f}  elongation={np.std(trajs[box_idx][:,0])/np.std(trajs[box_idx][:,1]):.1f}")
    print(f"  Loop: IC #{loop_idx}  R0={info[loop_idx][1]:.1f}  "
          f"f_tan={info[loop_idx][2]:.2f}  Lz_zc={np.sum(np.diff(np.sign(trajs[loop_idx][:,0]*trajs[loop_idx][:,4]-trajs[loop_idx][:,1]*trajs[loop_idx][:,3]))!=0)}")

    # ----- plot -----
    for idx, label, colour, fname in [
        (box_idx, 'Box (along bar, corotating)', 'blue', 'bar_box_orbit.pdf'),
        (loop_idx, 'Loop (corotating)', 'red', 'bar_loop_orbit.pdf')]:
        tr = trajs[idx]
        lim = np.max(np.abs(tr[:, :3]))
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
        for ax, (i, j), xl, yl in [
            (axes[0], (0, 1), 'x [kpc]', 'y [kpc]'),
            (axes[1], (0, 2), 'x [kpc]', 'z [kpc]'),
            (axes[2], (1, 2), 'y [kpc]', 'z [kpc]')]:
            ax.plot(tr[:, i], tr[:, j], '-', lw=0.5, color=colour)
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
            ax.set_aspect('equal'); ax.set_xlabel(xl); ax.set_ylabel(yl)
            # Mark bar boundary on x-y panel
            if ax == axes[0]:
                ax.add_patch(Ellipse((0, 0),
                    width=2*BAR_AXIS_X, height=2*BAR_AXIS_X*BAR_AXIS_RATIO,
                    edgecolor='grey', facecolor='none', ls='--', lw=1, alpha=0.6))
        axes[0].set_title(f'Barred (corotating): {label}')
        fig.tight_layout(); fig.savefig(fname); plt.close()
        print(f"  Saved {fname}")

    print("=== Demo 3 complete ===")

if __name__ == '__main__':
    main()
