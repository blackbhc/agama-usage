#!/usr/bin/env python3
"""
Demo 1: Exponential Disk + NFW Dark Matter Halo
Generate N-body initial conditions for a disk+halo galaxy model.
Reference: arXiv:2602.19995

Output: galaxy_ics.dat (columns: x y z vx vy vz mass)
"""
import agama, numpy as np, os

CONFIG = os.path.join(os.path.dirname(__file__), '..', 'configs',
                       'test1_exponential_disk_nfw.ini')
N_DISK = 20000
N_HALO = 5000
OUTPUT  = 'galaxy_ics.dat'

agama.setUnits(length=1, velocity=1, mass=1)

def main():
    print("=== Demo 1: Exp Disk + NFW DM Halo ICs ===")

    # 1. Load composite potential from INI
    pot = agama.Potential(file=CONFIG)
    vc = (-pot.force([[8, 0, 0]])[0, 0] * 8) ** 0.5
    print(f"  Host Vc(R=8) = {vc:.1f} km/s")

    # 2. Disk DF + sampling  (use composite potential)
    print(f"  Sampling {N_DISK} disk particles ...")
    df_disk = agama.DistributionFunction(
        type='quasiisothermal', potential=pot,
        mass=6e10, rdisk=3.0, hdisk=0.3,
        sigmar0=30.0, rsigmar=6.0, sigmamin=5.0)
    xv_d, m_d = agama.GalaxyModel(pot, df_disk).sample(N_DISK)
    print(f"    Disk: {len(xv_d)} particles, M={np.sum(m_d):.3e} Msun")

    # 3. Halo DF + sampling (use pure NFW for the Eddington inversion)
    #    The NFW mass=1.2e12 gives ~2.4e14 Msun enclosed within ~80 kpc
    print(f"  Sampling {N_HALO} halo particles ...")
    halo_pot = agama.Potential(type='NFW', mass=1.2e12, scaleRadius=15.0)
    df_halo = agama.DistributionFunction(type='quasispherical', potential=halo_pot)
    xv_h, m_h = agama.GalaxyModel(halo_pot, df_halo).sample(N_HALO)
    print(f"    Halo: {len(xv_h)} particles, M={np.sum(m_h):.3e} Msun")

    # 4. Combine and save as plain-text .dat
    xv = np.vstack([xv_d, xv_h])
    m  = np.concatenate([m_d, m_h])
    data = np.column_stack([xv, m])
    header = (f"x[kpc] y[kpc] z[kpc] vx[km/s] vy[km/s] vz[km/s] mass[Msun]\n"
              f"N={len(data)}  M_total={np.sum(m):.3e} Msun")
    np.savetxt(OUTPUT, data, fmt=['%12.6f']*6+['%14.6e'], header=header)
    print(f"  Saved: {OUTPUT}  ({len(data)} particles)")

if __name__ == '__main__':
    main()
