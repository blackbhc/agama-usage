#!/usr/bin/env python3
"""
Demo 4: Elliptical galaxy — box orbits and loop (tube) orbits
Triaxial Dehnen potential representing an elliptical galaxy.
Test various ICs, select one box orbit and one loop orbit.

Output: ellip_box_orbit.pdf, ellip_loop_orbit.pdf
"""
import agama, numpy as np, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'outputs')

agama.setUnits(length=1, velocity=1, mass=1)
TU = 0.978

def main():
    print("=== Demo 4: Elliptical galaxy — box & loop orbits ===")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Triaxial elliptical galaxy: Dehnen with axis ratios (1:0.8:0.6)
    pot = agama.Potential(
        type='Dehnen', mass=1e11, scaleRadius=2.0, gamma=1.0,
        axisRatioY=0.8, axisRatioZ=0.6)
    print("  Potential: triaxial Dehnen (1:0.8:0.6), M=1e11, a=2 kpc")

    # Generate many ICs at various positions and velocity directions
    np.random.seed(456)
    candidates = []
    for R in [1.5, 2.0, 3.0, 4.0, 5.0]:
        vc = (-R * pot.force([[R, 0, 0]])[0, 0])**0.5
        for e_factor in [0.3, 0.5, 0.7, 0.9]:
            v_phi = vc * e_factor
            for phi in [0, np.pi/3, 2*np.pi/3]:
                for tilt in [0, np.pi/6]:
                    x = R * np.cos(phi)
                    y = R * np.sin(phi) * 0.8
                    z = R * np.sin(tilt) * 0.3
                    vx = -v_phi * np.sin(phi)
                    vy =  v_phi * np.cos(phi) * 0.8
                    vz = v_phi * np.sin(tilt) * 0.2
                    candidates.append([x, y, z, vx, vy, vz])
    candidates = np.array(candidates)
    print(f"  Testing {len(candidates)} ICs ...")

    # Integrate
    T_int = 8.0 / TU
    _, trajs = agama.orbit(potential=pot, ic=candidates, time=T_int,
                           trajsize=2001, separateTime=True)

    # Classify by angular momentum behaviour
    box_idx, loop_idx = None, None
    for i in range(len(candidates)):
        tr = trajs[i]
        L = np.cross(tr[:, :3], tr[:, 3:6])
        Lz = L[:, 2]
        L_mag = np.mean(np.sqrt(np.sum(L**2, axis=1)))
        zc = np.sum(np.diff(np.sign(Lz)) != 0)
        # Box: signs of all 3 L components flip frequently
        # Loop (tube): at least one component preserves sign
        Lx_flip = np.sum(np.diff(np.sign(L[:, 0])) != 0)
        Ly_flip = np.sum(np.diff(np.sign(L[:, 1])) != 0)
        Lz_flip = zc
        all_flip = Lx_flip + Ly_flip + Lz_flip
        if all_flip > 60 and box_idx is None:
            box_idx = i
            print(f"  Box candidate {i}: flips={all_flip}")
        if all_flip < 20 and L_mag > 0.1 and loop_idx is None:
            loop_idx = i
            print(f"  Loop candidate {i}: flips={all_flip}")

    if box_idx is None or loop_idx is None:
        print("  Using fallback indices")
        box_idx = box_idx or 3
        loop_idx = loop_idx or 12

    # Plot box orbit
    tr = trajs[box_idx]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    axes[0].plot(tr[:, 0], tr[:, 1], 'b-', lw=0.5); axes[0].set_aspect('equal')
    axes[0].set_xlabel('x [kpc]'); axes[0].set_ylabel('y [kpc]')
    axes[0].set_title(f'Triaxial: Box Orbit (IC #{box_idx})')
    axes[1].plot(tr[:, 0], tr[:, 2], 'b-', lw=0.5); axes[1].set_aspect('equal')
    axes[1].set_xlabel('x [kpc]'); axes[1].set_ylabel('z [kpc]')
    axes[2].plot(tr[:, 1], tr[:, 2], 'b-', lw=0.5); axes[2].set_aspect('equal')
    axes[2].set_xlabel('y [kpc]'); axes[2].set_ylabel('z [kpc]')
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'ellip_box_orbit.pdf'))
    plt.close()
    print(f"  Saved ellip_box_orbit.pdf (IC #{box_idx})")

    # Plot loop (tube) orbit
    tr = trajs[loop_idx]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    axes[0].plot(tr[:, 0], tr[:, 1], 'r-', lw=0.5); axes[0].set_aspect('equal')
    axes[0].set_xlabel('x [kpc]'); axes[0].set_ylabel('y [kpc]')
    axes[0].set_title(f'Triaxial: Loop/Tube Orbit (IC #{loop_idx})')
    axes[1].plot(tr[:, 0], tr[:, 2], 'r-', lw=0.5); axes[1].set_aspect('equal')
    axes[1].set_xlabel('x [kpc]'); axes[1].set_ylabel('z [kpc]')
    axes[2].plot(tr[:, 1], tr[:, 2], 'r-', lw=0.5); axes[2].set_aspect('equal')
    axes[2].set_xlabel('y [kpc]'); axes[2].set_ylabel('z [kpc]')
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'ellip_loop_orbit.pdf'))
    plt.close()
    print(f"  Saved ellip_loop_orbit.pdf (IC #{loop_idx})")

    print("=== Demo 4 complete ===")

if __name__ == '__main__':
    main()
