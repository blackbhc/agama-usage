#!/usr/bin/env python3
"""
Demo 3: Barred potential -- box orbits and loop orbits in the corotating frame.
The Ferrers bar (axisRatioY=0.35) is elongated along x and rotates at Omega=-40.
In the corotating frame, the bar is stationary.
True box orbits oscillate along the bar (x) filling a rectangular region.
Loop orbits circulate around the center with a preserved Lz sign.

Output: bar_box_orbit.pdf, bar_loop_orbit.pdf (identical axis limits)
"""
import agama, numpy as np, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CONFIG = os.path.join(os.path.dirname(__file__), 'bar_potential.ini')
agama.setUnits(length=1, velocity=1, mass=1)
TU = 0.978

# Bar pattern speed in internal units (rad / time-unit)
Omega_bar = -40.0  # km/s/kpc = rad per internal time-unit

def main():
    print("=== Demo 3: Barred potential -- box & loop orbits (corotating frame) ===")
    print(f"  Omega_bar = {Omega_bar} rad/t-unit = {Omega_bar:.1f} km/s/kpc")

    # Read base potential from INI, then add rotation to the bar component
    # so that the bar rotates at Omega_bar in the inertial frame.
    base = agama.Potential(file=CONFIG)
    # The bar is component index 2 (0=bulge, 1=disk, 2=bar, 3=halo)
    # Add rotation modifier: bar rotates at Omega_bar
    pot_bar = agama.Potential(
        type='Ferrers', mass=2e10, scaleRadius=4.0,
        axisRatioY=0.35, axisRatioZ=0.2,
        rotation=[[0, 0], [1000, Omega_bar * 1000]])
    # Total potential: static components + rotating bar
    pot = agama.Potential(base[0], base[1], pot_bar, base[3])
    print("  Bar has rotation modifier (rotating in inertial frame)")

    # Generate ICs designed to produce box vs loop orbits in corotating frame
    # Box orbit: launch from bar major axis with very little tangential v
    # Loop orbit: enough tangential velocity to circulate
    candidates = []
    info = []

    # Box orbit candidates: from bar axis with little or no v_y
    for R in [2.5, 3.0, 3.5, 4.0, 4.5, 5.0]:
        vc = (-R * pot.force([[R, 0, 0]])[0, 0])**0.5
        # Very little tangential motion
        for f_tan in [0.001, 0.003, 0.005, 0.01, 0.03]:
            for f_rad in [0.01, 0.02, 0.05, 0.10, 0.20]:
                v_phi = vc * f_tan
                v_r = vc * f_rad * (1 if np.random.rand() > 0.5 else -1)
                candidates.append([R, 0, 0, v_r, v_phi, 0])
                info.append(('box', R, f_tan, f_rad))

    # Loop orbit candidates: from various positions with enough rotation
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
    print(f"  Testing {len(candidates)} ICs ...")

    # Integrate in corotating frame (Omega = bar pattern speed)
    T_int = 8.0 / TU  # ~8 Gyr
    _, trajs = agama.orbit(potential=pot, ic=candidates, time=T_int,
                           Omega=Omega_bar, trajsize=5001, separateTime=True)
    print("  Classification:")
    print("  | #  | type | R0 | f_tan | Lz_zc | x_std/y_std | class |")

    box_score = []; loop_score = []
    for i in range(len(candidates)):
        tr = trajs[i]
        x, y = tr[:, 0], tr[:, 1]
        Lz = x * tr[:, 4] - y * tr[:, 3]
        zc = np.sum(np.diff(np.sign(Lz)) != 0)
        x_std, y_std = np.std(x), np.std(y)
        aspect = min(x_std, y_std) / max(x_std, y_std) if max(x_std, y_std) > 0 else 1
        mean_Lz = np.mean(np.abs(Lz))

        # Box: many zero-crossings, fills both sides of bar (x-std > y-std for bar-aligned)
        is_box_like = zc > 30 and aspect > 0.3
        bs = zc * aspect if is_box_like else 0
        # Loop: few zero-crossings, significant circular motion
        is_loop_like = zc < 10 and mean_Lz > 10
        ls = -zc + mean_Lz * 0.1 if is_loop_like else 0

        box_score.append(bs)
        loop_score.append(ls)

        if i < 5 or is_box_like or is_loop_like:
            t = info[i][0]
            print(f"  | {i:3d} | {t:4s} | {info[i][1]:.1f} | {info[i][2]:.2f} | {zc:5d} | {aspect:.2f} | {'BOX->' if is_box_like else 'LOOP->' if is_loop_like else '     '} |")

    box_idx = np.argmax(box_score)
    loop_idx = np.argmax(loop_score)
    print(f"\n  Selected box:  IC #{box_idx} (R0={info[box_idx][1]:.1f}, f_tan={info[box_idx][2]:.2f})")
    print(f"  Selected loop: IC #{loop_idx} (R0={info[loop_idx][1]:.1f}, f_tan={info[loop_idx][2]:.2f})")

    # Compute Lz time series for box orbit to verify
    tr_b = trajs[box_idx]
    Lz_b = tr_b[:,0]*tr_b[:,4] - tr_b[:,1]*tr_b[:,3]
    print(f"  Box Lz: mean={np.mean(Lz_b):.1f}, std={np.std(Lz_b):.1f}, sign_changes={np.sum(np.diff(np.sign(Lz_b))!=0)}")

    # Plot limits per orbit
    lim_box = np.max(np.abs(trajs[box_idx][:, :3]))
    lim_loop = np.max(np.abs(trajs[loop_idx][:, :3]))

    for idx, label, colour, fname, lim in [
        (box_idx, 'Box (along bar, corotating)', 'blue', 'bar_box_orbit.pdf', lim_box),
        (loop_idx, 'Loop (corotating)', 'red', 'bar_loop_orbit.pdf', lim_loop)]:
        tr = trajs[idx]
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
        for ax, (i, j), xl, yl in [
            (axes[0], (0, 1), 'x [kpc]', 'y [kpc]'),
            (axes[1], (0, 2), 'x [kpc]', 'z [kpc]'),
            (axes[2], (1, 2), 'y [kpc]', 'z [kpc]')]:
            ax.plot(tr[:, i], tr[:, j], '-', lw=0.5, color=colour)
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
            ax.set_aspect('equal'); ax.set_xlabel(xl); ax.set_ylabel(yl)
        axes[0].set_title(f'Barred: {label}')
        fig.tight_layout(); fig.savefig(fname); plt.close()
        print(f"  Saved {fname}")

    print("=== Demo 3 complete ===")

if __name__ == '__main__':
    main()
