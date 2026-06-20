#!/Users/chenbinhui/miniconda3/bin/python
"""
Test 1: Exponential Disk + NFW Dark Matter Halo
Construct a realistic galaxy model with disk + dark halo.
Reference: arxiv:2602.19995
"""
import agama
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# === Step 1: Setup ===
agama.setUnits(length=1, velocity=1, mass=1)  # kpc, km/s, Msun
TIME_UNIT = 0.978  # Gyr per time unit

results_dir = os.path.dirname(os.path.abspath(__file__)) + '/results'
os.makedirs(results_dir, exist_ok=True)

print("=== Test 1: Exponential Disk + NFW Halo ===")
print(f"G = {agama.G:.6e}")

# === Step 2: Build Potentials ===
# Disk: exponential, M=5e10 Msun, R_d=3 kpc, h_z=0.3 kpc (positive = exponential)
disk_pot = {'type': 'Disk', 'mass': 5e10, 'scaleRadius': 3.0, 'scaleHeight': 0.3}

# Halo: NFW, M_200 ~ 1.2e12 Msun, R_s = 15 kpc
halo_pot = agama.Potential({'type': 'NFW', 'mass': 1.2e12, 'scaleRadius': 15.0})

# Composite potential
pot = agama.Potential(disk_pot, halo_pot)
print("Potentials created successfully")

# === Step 3: Rotation Curve ===
R_vals = np.logspace(-1, 1.5, 100)  # 0.1 to ~31.6 kpc
xyz = np.column_stack([R_vals, np.zeros_like(R_vals), np.zeros_like(R_vals)])

# Total force — pass the (N,3) array directly
force_total = pot.force(xyz)
Fx, Fy, Fz = force_total[:, 0], force_total[:, 1], force_total[:, 2]
vc_total = np.sqrt(-R_vals * Fx)

# Component decomposition - force from each component
force_disk = pot[0].force(xyz)
force_halo = pot[1].force(xyz)
Fx_disk, Fy_disk, Fz_disk = force_disk[:, 0], force_disk[:, 1], force_disk[:, 2]
Fx_halo, Fy_halo, Fz_halo = force_halo[:, 0], force_halo[:, 1], force_halo[:, 2]
vc_disk = np.sqrt(-R_vals * Fx_disk) * np.sign(Fx_disk)  # sign to handle negative
vc_halo = np.sqrt(-R_vals * Fx_halo) * np.sign(Fx_halo)

# Only show real values
vc_disk = np.where(np.isfinite(vc_disk) & (vc_disk > 0), vc_disk, 0)
vc_halo = np.where(np.isfinite(vc_halo) & (vc_halo > 0), vc_halo, 0)

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(R_vals, vc_total, 'k-', linewidth=2, label='Total')
ax.plot(R_vals, vc_disk, 'b--', label='Disk')
ax.plot(R_vals, vc_halo, 'r--', label='Halo')
ax.set_xlabel('R [kpc]')
ax.set_ylabel('V_c [km/s]')
ax.set_xscale('log')
ax.set_title('Rotation Curve: Exponential Disk + NFW Halo')
ax.legend()
ax.grid(True, alpha=0.3)
fig.savefig(f'{results_dir}/rotation_curve.pdf', bbox_inches='tight')
plt.close()
print(f"Rotation curve saved. V_c max ~ {np.max(vc_total):.0f} km/s at R ~ {R_vals[np.argmax(vc_total)]:.1f} kpc")

# === Step 4: Create Distribution Functions ===
# Halo DF: use the spherical NFW-only potential for Eddington inversion
# then combine with composite potential in GalaxyModel
df_halo = agama.DistributionFunction(type='quasispherical',
    potential=halo_pot)
print("Halo DF created (QuasiSpherical, isotropic)")

# Disk DF: QuasiIsothermal
R_d = 3.0
Sigma0 = 5e10 / (2 * np.pi * R_d**2)  # for density profile reference
df_disk = agama.DistributionFunction(type='quasiisothermal',
    potential=pot,
    mass=5e10, rdisk=R_d, hdisk=0.3,
    sigmar0=30.0, rsigmar=3.0, sigmamin=5.0)
print(f"Disk DF created (QuasiIsothermal), Sigma0={Sigma0:.2e}")

# === Step 5: Sample Particles ===
# Halo particles: sample from spherical potential, then use composite for consistency
gm_halo = agama.GalaxyModel(halo_pot, df_halo)
print("Sampling 10000 halo particles (from spherical NFW potential)...")
xv_halo, m_halo = gm_halo.sample(10000)
print(f"  Halo: {len(xv_halo)} particles, total mass = {np.sum(m_halo):.2e} Msun")

# Disk particles
gm_disk = agama.GalaxyModel(pot, df_disk)
print("Sampling 5000 disk particles...")
xv_disk, m_disk = gm_disk.sample(5000)
print(f"  Disk: {len(xv_disk)} particles, total mass = {np.sum(m_disk):.2e} Msun")

# Combine
xv_all = np.vstack([xv_halo, xv_disk])
m_all = np.concatenate([m_halo, m_disk])

# === Step 6: Density Profile Verification ===
# Compute radial distances
r_halo = np.sqrt(np.sum(xv_halo[:, :3]**2, axis=1))
r_disk = np.sqrt(np.sum(xv_disk[:, :3]**2, axis=1))
R_disk = np.sqrt(xv_disk[:, 0]**2 + xv_disk[:, 1]**2)

# Halo 3D density profile
r_bins = np.logspace(-0.5, 2, 50)
r_center = np.sqrt(r_bins[:-1] * r_bins[1:])
r_width = np.diff(r_bins)
rho_bins = np.zeros(len(r_bins) - 1)
for i in range(len(r_bins) - 1):
    mask = (r_halo >= r_bins[i]) & (r_halo < r_bins[i+1])
    if np.sum(mask) > 0:
        vol = 4/3 * np.pi * (r_bins[i+1]**3 - r_bins[i]**3)
        rho_bins[i] = np.sum(m_halo[mask]) / vol

# Disk surface density profile
R_bins = np.linspace(0, 20, 41)
R_center_sd = (R_bins[1:] + R_bins[:-1]) / 2
Sigma_bins = np.zeros(len(R_bins) - 1)
for i in range(len(R_bins) - 1):
    mask = (R_disk >= R_bins[i]) & (R_disk < R_bins[i+1])
    if np.sum(mask) > 0:
        area = np.pi * (R_bins[i+1]**2 - R_bins[i]**2)
        Sigma_bins[i] = np.sum(m_disk[mask]) / area

# Plot density profiles
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Halo density
mask_good = rho_bins > 0
ax1.loglog(r_center[mask_good], rho_bins[mask_good], 'bo', markersize=4, label='Sampled')
# NFW reference
R_s = 15.0
M_200 = 1.2e12
rho_0 = M_200 / (4 * np.pi * R_s**3 * (np.log(1+1) - 1/(1+1)))
r_model = r_center
rho_nfw = rho_0 / ((r_model/R_s) * (1 + r_model/R_s)**2)
ax1.loglog(r_model, rho_nfw, 'r-', label='NFW (input)')
ax1.set_xlabel('r [kpc]')
ax1.set_ylabel(r'$\rho$ [M$_\odot$/kpc$^3$]')
ax1.set_title('Halo Density Profile')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Disk surface density
mask_good_sd = Sigma_bins > 0
ax2.semilogy(R_center_sd[mask_good_sd], Sigma_bins[mask_good_sd], 'bo', markersize=4, label='Sampled')
Sigma_model = Sigma0 * np.exp(-R_center_sd / R_d)
ax2.semilogy(R_center_sd, Sigma_model, 'r-', label=f'exp(-R/{R_d} kpc)')
ax2.set_xlabel('R [kpc]')
ax2.set_ylabel(r'$\Sigma$ [M$_\odot$/kpc$^2$]')
ax2.set_title('Disk Surface Density')
ax2.legend()
ax2.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(f'{results_dir}/density_profiles.pdf', bbox_inches='tight')
plt.close()
print("Density profile plots saved")

# === Step 7: Save ICs ===
# Combined set
agama.writeSnapshot(f'{results_dir}/galaxy_ICs.gadget', (xv_all, m_all), 'gadget')
# Component-split
agama.writeSnapshot(f'{results_dir}/halo_ICs.gadget', (xv_halo, m_halo), 'gadget')
agama.writeSnapshot(f'{results_dir}/disk_ICs.gadget', (xv_disk, m_disk), 'gadget')
print("ICs saved to gadget format")

# === Step 8: Quick Orbit Check ===
# Integrate a test particle at R=8 kpc with circular velocity
R0 = 8.0
vc_at_R0 = np.sqrt(-R0 * pot.force(np.array([[R0, 0, 0]]))[0, 0])
print(f"\nV_c(R=8 kpc) = {vc_at_R0:.1f} km/s")

# Expected circular period: T = 2*pi*R / V_c
T_exp = 2 * np.pi * R0 / vc_at_R0  # in time units
print(f"Expected circular period ~ {T_exp:.1f} (~ {T_exp * TIME_UNIT:.2f} Gyr)")

ic = [R0, 0, 0, 0, vc_at_R0, 0]
times, traj = agama.orbit(potential=pot, ic=ic, time=3*T_exp, trajsize=3001)

# Plot orbit
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(traj[:, 0], traj[:, 1], 'b-', linewidth=0.5)
ax1.plot(traj[0, 0], traj[0, 1], 'go', label='Start', markersize=8)
ax1.plot(traj[-1, 0], traj[-1, 1], 'ro', label='End', markersize=8)
ax1.set_xlabel('x [kpc]')
ax1.set_ylabel('y [kpc]')
ax1.set_title(f'Orbit at R={R0} kpc (V_c={vc_at_R0:.0f} km/s)')
ax1.set_aspect('equal')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Radial range
R_orbit = np.sqrt(traj[:, 0]**2 + traj[:, 1]**2)
ax2.plot(times * TIME_UNIT, R_orbit, 'k-', linewidth=0.8)
ax2.set_xlabel('Time [Gyr]')
ax2.set_ylabel('R [kpc]')
ax2.set_title('Radial Distance vs Time')
ax2.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(f'{results_dir}/test_orbit.pdf', bbox_inches='tight')
plt.close()
print("Test orbit plot saved")

# Energy conservation check
E = pot.potential(traj[:, :3]) + \
    0.5 * (traj[:, 3]**2 + traj[:, 4]**2 + traj[:, 5]**2)
dE = (E - E[0]) / np.abs(E[0])
print(f"Energy conservation: max|dE/E| = {np.max(np.abs(dE)):.2e}")
assert np.max(np.abs(dE)) < 1e-7, "Energy conservation FAILED!"

# Also plot particle distribution
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.scatter(xv_halo[:, 0], xv_halo[:, 1], s=0.3, alpha=0.5, c='r', label='Halo')
ax1.scatter(xv_disk[:, 0], xv_disk[:, 1], s=0.3, alpha=0.5, c='b', label='Disk')
ax1.set_xlabel('x [kpc]')
ax1.set_ylabel('y [kpc]')
ax1.set_title('Particle Distribution (xy plane)')
ax1.set_aspect('equal')
ax1.legend(markerscale=5)
ax1.grid(True, alpha=0.3)

ax2.scatter(xv_halo[:, 0], xv_halo[:, 2], s=0.3, alpha=0.5, c='r', label='Halo')
ax2.scatter(xv_disk[:, 0], xv_disk[:, 2], s=0.3, alpha=0.5, c='b', label='Disk')
ax2.set_xlabel('x [kpc]')
ax2.set_ylabel('z [kpc]')
ax2.set_title('Particle Distribution (xz plane)')
ax2.legend(markerscale=5)
ax2.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(f'{results_dir}/particle_distribution.pdf', bbox_inches='tight')
plt.close()
print("Particle distribution plot saved")

print("\n=== Test 1 COMPLETE ===")
print(f"All outputs saved to {results_dir}/")