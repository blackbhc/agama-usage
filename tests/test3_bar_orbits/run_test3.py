#!/Users/chenbinhui/miniconda3/bin/python
"""
Test 3: Barred Potential Orbit Integration (20 particles, 10 Gyr)
20 test particles in a barred potential, 5 figures × 4 orbits each.
"""
import agama
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

agama.setUnits(length=1, velocity=1, mass=1)  # kpc, km/s, Msun
TIME_UNIT = 0.978  # Gyr per time unit

results_dir = os.path.dirname(os.path.abspath(__file__)) + '/results'
os.makedirs(results_dir, exist_ok=True)

print("=== Test 3: Barred Potential Orbit Integration ===")

# === Step 1: Create Barred Potential ===
# Axisymmetric components
bulge = {'type': 'Dehnen', 'mass': 2e10, 'scaleRadius': 0.5, 'gamma': 0}
disk = {'type': 'Disk', 'mass': 5e10, 'scaleRadius': 3.0, 'scaleHeight': 0.3}
halo = {'type': 'Logarithmic', 'v0': 200.0, 'scaleRadius': 5.0}

# Bar: Ferrers triaxial ellipsoid
bar = {'type': 'Ferrers', 'mass': 1.5e10, 'scaleRadius': 5.0,
       'axisRatioY': 0.3, 'axisRatioZ': 0.2}

pot = agama.Potential(bulge, disk, halo, bar)

# Pattern speed (km/s/kpc, negative = clockwise)
Omega_bar = -40.0
print(f"Bar pattern speed: {Omega_bar} km/s/kpc")
print(f"Corotation radius ~ {200.0/abs(Omega_bar):.1f} kpc (approx)")

# === Step 2: Generate 20 Initial Conditions ===
np.random.seed(42)  # reproducible
ics = []
for i in range(20):
    R = 2.0 + np.random.rand() * 6.0  # 2-8 kpc
    phi = np.random.rand() * 2 * np.pi
    vc = np.sqrt(-R * pot.force(np.array([[R, 0, 0]]))[0, 0])
    vr = (np.random.rand() - 0.5) * 40
    vt = vc + (np.random.rand() - 0.5) * 30
    x = R * np.cos(phi)
    y = R * np.sin(phi)
    vx = vr * np.cos(phi) - vt * np.sin(phi)
    vy = vr * np.sin(phi) + vt * np.cos(phi)
    ics.append([x, y, 0, vx, vy, 0])
ics = np.array(ics)
print(f"Generated {len(ics)} initial conditions")

# === Step 3: Integrate Orbits for 10 Gyr ===
t_Gyr = 10.0
t_internal = t_Gyr / TIME_UNIT  # convert to internal time units
print(f"Integrating for {t_Gyr} Gyr ({t_internal:.1f} internal units)")

times_all, trajs_all = agama.orbit(
    potential=pot, ic=ics, time=t_internal,
    Omega=Omega_bar, trajsize=10000, separateTime=True)
# trajs_all shape: (20, 10000, 6)

print(f"Orbits integrated: times shape {np.shape(times_all)}, "
      f"trajs shape {np.shape(trajs_all)}")

# === Step 4: Plot 5 Figures × 4 Particles Each ===
for fig_idx in range(5):
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    for sub_idx, ax in enumerate(axes.flat):
        p_idx = fig_idx * 4 + sub_idx
        if p_idx >= 20:
            ax.text(0.5, 0.5, 'N/A', transform=ax.transAxes, ha='center')
            continue
        traj = trajs_all[p_idx]
        ax.plot(traj[:, 0], traj[:, 1], 'b-', linewidth=0.3, alpha=0.8)
        ax.plot(ics[p_idx, 0], ics[p_idx, 1], 'go', markersize=6, label='Start')
        ax.plot(traj[-1, 0], traj[-1, 1], 'ro', markersize=6, label='End')
        ax.set_xlim(-10, 10)
        ax.set_ylim(-10, 10)
        ax.set_aspect('equal')
        ax.set_xlabel('x [kpc]')
        ax.set_ylabel('y [kpc]')
        ax.set_title(f'Particle {p_idx+1}: R₀={np.sqrt(ics[p_idx,0]**2+ics[p_idx,1]**2):.1f} kpc')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    fig.suptitle(f'Orbits in Barred Potential (Omega={Omega_bar} km/s/kpc), Figure {fig_idx+1}/5',
                 fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(f'{results_dir}/bar_orbits_fig{fig_idx+1}.pdf', bbox_inches='tight')
    plt.close()

print(f"5 figures saved to {results_dir}/")

# === Step 5: Quick Analysis — Energy conservation ===
# In rotating frame, Jacobi integral is conserved
traj_p0 = trajs_all[0]  # (10000, 6)
E_jacobi = pot.potential(traj_p0[:, :3]) + \
    0.5 * (traj_p0[:, 3]**2 + traj_p0[:, 4]**2 + traj_p0[:, 5]**2) \
    - 0.5 * (Omega_bar * np.sqrt(traj_p0[:, 0]**2 + traj_p0[:, 1]**2))**2
dE = (E_jacobi - E_jacobi[0]) / np.abs(E_jacobi[0])
print(f"Particle 1 Jacobi integral conservation: max|dE/E| = {np.max(np.abs(dE)):.2e}")

print("\n=== Test 3 COMPLETE ===")
print(f"All outputs saved to {results_dir}/")