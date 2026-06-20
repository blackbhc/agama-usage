#!/usr/bin/env python3
"""
Demo 1: Exponential Disk + NFW Dark Matter Halo
Generate N-body initial conditions for a disk+halo galaxy model.
Reference: arXiv:2602.19995

Outputs:
  disk_particles.dat  — 20000 disk particles (x y z vx vy vz mass)
  halo_particles.dat  — 5000  halo particles (x y z vx vy vz mass)
"""
import agama, numpy as np, os

CONFIG = os.path.join(os.path.dirname(__file__), 'galaxy.ini')
N_DISK = 20000
N_HALO = 5000

agama.setUnits(length=1, velocity=1, mass=1)

def main():
    print("=== Demo 1: Exp Disk + NFW DM Halo ICs ===")

    # 1. Composite potential from INI
    pot = agama.Potential(file=CONFIG)
    vc = (-pot.force([[8, 0, 0]])[0, 0] * 8)**0.5
    print(f"  Host Vc(R=8) = {vc:.1f} km/s")

    # 2. Disk DF + sample (in the full composite potential)
    print(f"  Sampling {N_DISK} disk particles ...")
    df_disk = agama.DistributionFunction(
        type='quasiisothermal', potential=pot,
        mass=6e10, rdisk=3.0, hdisk=0.3,
        sigmar0=30.0, rsigmar=6.0, sigmamin=5.0)
    xv_d, m_d = agama.GalaxyModel(pot, df_disk).sample(N_DISK)
    print(f"    Disk: {len(xv_d)} particles, M={np.sum(m_d):.3e} Msun")
    data_d = np.column_stack([xv_d, m_d])
    np.savetxt('disk_particles.dat', data_d,
               fmt=['%12.6f']*6+['%14.6e'],
               header=f'x[kpc] y[kpc] z[kpc] vx[km/s] vy[km/s] vz[km/s] mass[Msun]')
    print(f"  Saved: disk_particles.dat")

    # 3. Halo: use pure NFW for the Eddington inversion (ref: example_self_consistent_model3.py)
    print(f"  Sampling {N_HALO} halo particles ...")
    halo_pot = agama.Potential(type='NFW', mass=1.2e12, scaleRadius=15.0)
    df_halo = agama.DistributionFunction(type='quasispherical', potential=halo_pot)
    xv_h, m_h = agama.GalaxyModel(halo_pot, df_halo).sample(N_HALO)
    print(f"    Halo: {len(xv_h)} particles, M={np.sum(m_h):.3e} Msun")
    data_h = np.column_stack([xv_h, m_h])
    np.savetxt('halo_particles.dat', data_h,
               fmt=['%12.6f']*6+['%14.6e'],
               header=f'x[kpc] y[kpc] z[kpc] vx[km/s] vy[km/s] vz[km/s] mass[Msun]')
    print(f"  Saved: halo_particles.dat")

    print("=== Demo 1 complete ===")

if __name__ == '__main__':
    main()
