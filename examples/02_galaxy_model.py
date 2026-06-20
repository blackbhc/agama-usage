#!/usr/bin/env python3
"""Example 2: Galaxy model — disk DF, sampling, moments, snapshot I/O."""
import agama, numpy as np, os

config_dir = os.path.join(os.path.dirname(__file__), '..', 'configs')
agama.setUnits(length=1, velocity=1, mass=1)

print("=" * 60)
print("1. Load potential from INI...")
pot = agama.Potential(file=os.path.join(config_dir, 'test1_exponential_disk_nfw.ini'))
vc = (-pot.force([[8, 0, 0]])[0, 0] * 8) ** 0.5
print(f"   Vc(R=8) = {vc:.1f} km/s")

print("2. Create QuasiIsothermal disk DF...")
df = agama.DistributionFunction(type='quasiisothermal', potential=pot,
    mass=5e10, rdisk=3.0, hdisk=0.3,
    sigmar0=30.0, rsigmar=6.0, sigmamin=5.0)
print("   Disk DF OK")

print("3. Sample 2000 disk particles...")
xv, m = agama.GalaxyModel(pot, df).sample(2000)
print(f"   {len(xv)} particles, mass={np.sum(m):.3e} Msun")

print("4. Velocity moments at R=[4,8,12] kpc...")
xyz = np.column_stack([[4.0, 8.0, 12.0], np.zeros(3), np.zeros(3)])
rho, vel, vel2 = agama.GalaxyModel(pot, df).moments(xyz, dens=True, vel=True, vel2=True)
for i, R in enumerate([4.0, 8.0, 12.0]):
    sr = np.sqrt(vel2[i, 0])   # sigma_xx
    sz = np.sqrt(vel2[i, 2])   # sigma_zz
    print(f"   R={R:5.1f}: v_phi={vel[i,1]:5.0f}, sr={sr:4.0f}, sz={sz:4.0f} km/s")

print("5. Save/load snapshot...")
snap = os.path.join(os.path.dirname(__file__), '..', 'galaxy_ICs.gadget')
agama.writeSnapshot(snap, (xv, m), 'gadget')
posvel, mass = agama.readSnapshot(snap)
assert len(posvel) == len(xv)
print(f"   Verified: {len(posvel)} particles round-trip OK")

print("\n✅ Example 2 complete")
