#!/usr/bin/env python3
"""Example 3: Orbit integration in various potentials."""
import agama, numpy as np, os

config_dir = os.path.join(os.path.dirname(__file__), '..', 'configs')
agama.setUnits(length=1, velocity=1, mass=1)
TU = 0.978  # Gyr per time unit

print("=" * 60)
print("1. Simple orbit (axisymmetric MW-like)...")
pot = agama.Potential(file=os.path.join(config_dir, 'test2_mw_host.ini'))
vc = (-pot.force([[8, 0, 0]])[0, 0] * 8) ** 0.5
t, traj = agama.orbit(potential=pot, ic=[8, 0, 0, 0, vc, 15],
                       time=5/TU, trajsize=501, separateTime=True)
R = np.sqrt(traj[:, 0]**2 + traj[:, 1]**2)
E = pot.potential(traj[:, :3]) + 0.5 * np.sum(traj[:, 3:]**2, axis=1)
dE = np.max(np.abs((E - E[0]) / E[0]))
print(f"   R range={np.min(R):.2f}-{np.max(R):.2f} kpc, dE/E={dE:.2e}")
assert dE < 1e-7

print("2. Rotating barred potential...")
pot_b = agama.Potential(file=os.path.join(config_dir, 'test3_bar_potential.ini'))
t, traj = agama.orbit(potential=pot_b, ic=[5, 0, 0, 0, 180, 0],
                       Omega=-40, time=3/TU, trajsize=2001, separateTime=True)
print(f"   x range=[{np.min(traj[:,0]):.1f}, {np.max(traj[:,0]):.1f}]")

print("3. Batch integration (10 orbits)...")
R_v = np.random.uniform(3, 12, 10)
vcs = np.array([(-R * pot.force([[R, 0, 0]])[0, 0])**0.5 for R in R_v])
ics = np.column_stack([R_v, np.zeros(10), np.zeros(10), np.zeros(10), vcs, np.zeros(10)])
t_arr, tr_arr = agama.orbit(potential=pot, ic=ics, time=2/TU, trajsize=201, separateTime=True)
print(f"   Shape: {tr_arr.shape}")

pot2 = agama.Potential(file=os.path.join(config_dir, "test1_exponential_disk_nfw.ini"))
print("4. Actions and frequencies...")
af = agama.ActionFinder(pot2)
J = af(ics)
print(f"   Jr: [{np.min(J[:,0]):.0f}, {np.max(J[:,0]):.0f}] kpc*km/s")

print("5. Time-dependent (moving satellite)...")
sat = agama.Potential(type='Plummer', mass=1e9, scaleRadius=0.5,
                       center=[[0, 50, 0, 0], [3/TU, 0, 0, 0]])
pot_td = agama.Potential(pot, sat)
t, traj = agama.orbit(potential=pot_td, ic=[45, 0, 0, 0, 180, 0],
                       time=2/TU, trajsize=201, separateTime=True)
R = np.sqrt(traj[:, 0]**2 + traj[:, 1]**2)
print(f"   R range=[{np.min(R):.1f}, {np.max(R):.1f}] kpc")

print("\n✅ Example 3 complete")
