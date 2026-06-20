# Test 3: Barred Potential Orbit Integration (20 particles, 10 Gyr)

**Scenario**: 20 test particles in a barred spiral potential, integrated for 10 Gyr in the co-rotating frame, plotted as 5 figures × 4 particles each.

## Steps
1. **Host potential**: bulge (Dehnen) + disk (Disk) + halo (Logarithmic)
2. **Bar**: Ferrers triaxial ellipsoid with pattern speed Omega=-40 km/s/kpc
3. **20 initial conditions**: randomly distributed in R=[3,10] kpc, with small velocity perturbations
4. **Integrate**: 10 Gyr with Omega=Omega_bar in co-rotating frame, 10000 output points
5. **Plot**: 5 figures, each with 4 particle orbits in xy-plane