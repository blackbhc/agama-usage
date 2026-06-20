#!/Users/chenbinhui/miniconda3/bin/python
"""
Test 2: 1000-Particle Satellite Merger with MW-like Galaxy
Simulate a satellite galaxy on an eccentric orbit around a Milky Way-like host.
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

print("=== Test 2: 1000-Particle Satellite Merger ===")

# === Step 1: MW Host Potential ===
# Bulge
bulge = {'type': 'Dehnen', 'mass': 2e10, 'scaleRadius': 0.5, 'gamma': 0}
# Disk
disk = {'type': 'Disk', 'mass': 5e10, 'scaleRadius': 3.0, 'scaleHeight': 0.3}
# Halo (Logarithmic for flat rotation curve)
halo = {'type': 'Logarithmic', 'v0': 200.0, 'scaleRadius': 5.0}

pot_host = agama.Potential(bulge, disk, halo)
print("Host potential created")

# Rotation curve check
R_vals = np.linspace(0.1, 50, 100)
xyz = np.column_stack([R_vals, np.zeros_like(R_vals), np.zeros_like(R_vals)])
force = pot_host.force(xyz)
vc = np.sqrt(-R_vals * force[:, 0])
print(f"V_c(R=8) = {np.interp(8.0, R_vals, vc):.1f} km/s, "
      f"V_c(R=50) = {np.interp(50.0, R_vals, vc):.1f} km/s")

# === Step 2: Satellite ===
R0_sat = 50.0  # kpc, initial position
Vcirc_sat = np.interp(R0_sat, R_vals, vc)
V0_sat = 0.7 * Vcirc_sat  # 70% circular = eccentric orbit

# Satellite potential (Plummer)
pot_sat = agama.Potential(type='Plummer', mass=1e9, scaleRadius=0.5)

# Satellite DF + sample
df_sat = agama.DistributionFunction(type='quasispherical', potential=pot_sat)
gm_sat = agama.GalaxyModel(pot_sat, df_sat)
xv_sat, m_sat = gm_sat.sample(1000)
print(f"Satellite: {len(xv_sat)} particles, mass sum = {np.sum(m_sat):.2e} Msun")

# Place satellite at its orbital position
xv_sat[:, 0] += R0_sat
xv_sat[:, 3] += V0_sat

# === Step 3: Integrate satellite center of mass ===
ic_sat = [R0_sat, 0, 0, 0, V0_sat, 0]
T_orb = 2 * np.pi * R0_sat / Vcirc_sat
t_total = 4.0 * T_orb  # 4 orbital periods

times, traj_sat = agama.orbit(potential=pot_host, ic=ic_sat,
                               time=t_total, trajsize=500)
print(f"Orbit integrated: {t_total:.1f} time units ({t_total*TIME_UNIT:.1f} Gyr)")
print(f"  Period ~ {T_orb:.1f} ({T_orb*TIME_UNIT:.1f} Gyr)")

# Plot orbit
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
ax1.plot(traj_sat[:, 0], traj_sat[:, 1], 'k-', linewidth=1)
ax1.plot(traj_sat[0, 0], traj_sat[0, 1], 'go', label='Start', markersize=8)
ax1.set_xlabel('x [kpc]'); ax1.set_ylabel('y [kpc]')
ax1.set_title('Satellite Center-of-Mass Orbit')
ax1.set_aspect('equal'); ax1.legend(); ax1.grid(True, alpha=0.3)

R_orbit = np.sqrt(traj_sat[:, 0]**2 + traj_sat[:, 1]**2)
ax2.plot(times * TIME_UNIT, R_orbit, 'k-')
ax2.axhline(R0_sat, color='g', linestyle='--', alpha=0.5, label='R_initial')
ax2.set_xlabel('Time [Gyr]'); ax2.set_ylabel('R [kpc]')
ax2.set_title('Orbital Radius vs Time'); ax2.legend(); ax2.grid(True, alpha=0.3)
fig.tight_layout(); fig.savefig(f'{results_dir}/satellite_orbit.pdf', bbox_inches='tight')
plt.close()

# === Step 4: Integrate test particles ===
# For each timestep, offset satellite particles to follow the COM trajectory
# Then integrate all particles in the host potential
n_particles = len(xv_sat)
n_snapshots = min(100, len(times))

# Sample timesteps
snap_indices = np.linspace(0, len(times)-1, n_snapshots, dtype=int)

# Store particle positions at selected snapshots
particle_snapshots = np.zeros((n_snapshots, n_particles, 6))
particle_snapshots[0] = xv_sat

# Integrate all particles in the static host potential
# (test particle approximation - no self-gravity of satellite)
t_particle = t_total
orbit_result = agama.orbit(potential=pot_host, ic=xv_sat,
                           time=t_particle, trajsize=n_snapshots,
                           separateTime=True)
# With separateTime=True: (N_times,) + (N_particles, N_times, 6)
times_p = orbit_result[0]  # (N_times,) array or (N_particles, N_times)
trajs_p = orbit_result[1]  # (N_particles, N_times, 6) with separateTime=True
# times_p is (N_particles, N_times) - use first particle's times
times_p_1d = times_p[0]

print(f"Integrated {n_particles} particles, {n_snapshots} snapshots")

# === Step 5: Analyze tidal stripping ===
# For each snapshot, compute the distance from satellite center
# trajs_p index: [particle, time, coord]
R_snap = np.sqrt(traj_sat[snap_indices, 0]**2 + traj_sat[snap_indices, 1]**2)

# Count bound particles (within tidal radius ~ R_satellite)
R_tidal = 2.0  # kpc, approximate tidal radius
bound_fraction = np.zeros(n_snapshots)
for i in range(n_snapshots):
    dx = trajs_p[:, i, 0] - traj_sat[snap_indices[i], 0]
    dy = trajs_p[:, i, 1] - traj_sat[snap_indices[i], 1]
    dz = trajs_p[:, i, 2] - traj_sat[snap_indices[i], 2]
    dr = np.sqrt(dx**2 + dy**2 + dz**2)
    bound_fraction[i] = np.sum(dr < R_tidal) / n_particles

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(times[snap_indices] * TIME_UNIT, bound_fraction * 100, 'b-', linewidth=1.5)
ax.set_xlabel('Time [Gyr]'); ax.set_ylabel('Bound Fraction [%]')
ax.set_title('Satellite Mass Loss (tidal stripping)')
ax.grid(True, alpha=0.3)
fig.savefig(f'{results_dir}/mass_loss.pdf', bbox_inches='tight')
plt.close()
print(f"Final bound fraction: {bound_fraction[-1]*100:.1f}%")

# === Step 6: Visualize particle distribution at snapshots ===
snap_to_plot = [0, n_snapshots//4, n_snapshots//2, -1]
labels = [f't={times[snap_indices[i]]*TIME_UNIT:.1f} Gyr' for i in snap_to_plot]

fig, axes = plt.subplots(2, 2, figsize=(12, 12))
for idx, ax in enumerate(axes.flat):
    i = snap_to_plot[idx]
    ax.scatter(trajs_p[:, i, 0], trajs_p[:, i, 1], s=1, alpha=0.5, c='b')
    ax.plot(traj_sat[snap_indices[i], 0], traj_sat[snap_indices[i], 1],
            'ro', markersize=5, label='Sat center')
    ax.set_xlabel('x [kpc]'); ax.set_ylabel('y [kpc]')
    ax.set_title(labels[idx])
    ax.set_aspect('equal'); ax.legend()
    ax.set_xlim(-60, 60); ax.set_ylim(-60, 60)

fig.tight_layout(); fig.savefig(f'{results_dir}/particle_snapshots.pdf', bbox_inches='tight')
plt.close()
print("Particle snapshot plot saved")

# Also plot final particle positions zoomed in
fig, ax = plt.subplots(figsize=(8, 8))
i_final = -1
ax.scatter(trajs_p[:, i_final, 0], trajs_p[:, i_final, 1], s=2, alpha=0.6, c='b')
ax.plot(traj_sat[snap_indices[i_final], 0], traj_sat[snap_indices[i_final], 1],
        'ro', markersize=8, label='Satellite center')
ax.set_xlabel('x [kpc]'); ax.set_ylabel('y [kpc]')
ax.set_title(f'Final Particle Distribution (t={times[-1]*TIME_UNIT:.1f} Gyr)')
ax.set_aspect('equal'); ax.legend(); ax.grid(True, alpha=0.3)
ax.set_xlim(-20, 20); ax.set_ylim(-20, 20)
fig.savefig(f'{results_dir}/final_distribution.pdf', bbox_inches='tight')
plt.close()

print("\n=== Test 2 COMPLETE ===")
print(f"All outputs saved to {results_dir}/")