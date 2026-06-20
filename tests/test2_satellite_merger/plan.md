# Test 2: 1000-Particle Satellite Merger with MW-like Galaxy

**Scenario**: A 1000-particle satellite orbiting and merging with a Milky Way-like galaxy.

## Steps
1. **MW Host Potential**: bulge (Spheroid) + disk (Disk) + halo (NFW/Logarithmic)
2. **Satellite**: Plummer sphere (mass=1e9 Msun, scaleRadius=0.5 kpc)
3. **Sample satellite particles**: 1000 from QuasiSpherical DF
4. **Place on eccentric orbit**: R_apocenter=50 kpc, 70% circular velocity
5. **Integrate**: orbit of satellite center + test particles in host potential
6. **Visualize**: particle snapshots at different times, orbital decay