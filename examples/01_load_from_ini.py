#!/usr/bin/env python3
"""Example 1: Load galaxy potentials from INI configuration files."""
import agama, numpy as np, os

config_dir = os.path.join(os.path.dirname(__file__), '..', 'configs')
agama.setUnits(length=1, velocity=1, mass=1)

print("=" * 60)
print("1. Loading composite potential from INI...")
pot = agama.Potential(file=os.path.join(config_dir, 'test2_mw_host.ini'))
print(f"   Components: {len(pot)}")
vc = (-pot.force([[8, 0, 0]])[0, 0] * 8) ** 0.5
print(f"   Vc(R=8) = {vc:.1f} km/s")

print("\n2. Loading McMillan17 MW model...")
pot2 = agama.Potential(file=os.path.join(config_dir, 'McMillan17.ini'))
print(f"   Components: {len(pot2)}")
vc2 = (-pot2.force([[8, 0, 0]])[0, 0] * 8) ** 0.5
print(f"   Vc(R=8) = {vc2:.1f} km/s")

print("\n3. Loading test1 (disk+NFW)...")
pot3 = agama.Potential(file=os.path.join(config_dir, 'test1_exponential_disk_nfw.ini'))
vc3 = (-pot3.force([[8, 0, 0]])[0, 0] * 8) ** 0.5
print(f"   Vc(R=8) = {vc3:.1f} km/s")

print("\n4. Loading Cautun20 (Multipole expansion)...")
pot4 = agama.Potential(file=os.path.join(config_dir, 'Cautun20.ini'))
vc4 = (-pot4.force([[8, 0, 0]])[0, 0] * 8) ** 0.5
print(f"   Vc(R=8) = {vc4:.1f} km/s")

print("\n5. Rotation curve (test2_mw_host)")
R = np.logspace(-1, 1.5, 20)
xyz = np.column_stack([R, np.zeros_like(R), np.zeros_like(R)])
fini = pot.force(xyz)
vc_ini = np.sqrt(-R * fini[:, 0])
for r, v in zip(R[::4], vc_ini[::4]):
    print(f"   R={r:6.2f} kpc, Vc={v:7.1f} km/s")

print("\n✅ Example 1 complete")
