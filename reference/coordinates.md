# Coordinate Systems Reference

AGAMA provides comprehensive coordinate transformations between Galactocentric, Galactic (l,b), and ICRS (RA,Dec) systems.

## 1. Galactocentric ↔ Galactic

```python
# Galactocentric → Galactic (l, b, distance, proper motion, vlos)
l, b, d, pml, pmb, vlos = agama.getGalacticFromGalactocentric(
    x, y, z, vx, vy, vz,
    galcen_distance=8.122,        # distance to GC [kpc]
    galcen_v_sun=(12.9, 245.6, 7.78),  # U, V, W solar motion [km/s]
    z_sun=0.0208)                 # solar height [kpc]

# Galactic → Galactocentric
pos = agama.getGalactocentricFromGalactic(
    l, b, d, pml, pmb, vlos,
    galcen_distance=8.122,
    galcen_v_sun=(12.9, 245.6, 7.78),
    z_sun=0.0208)
# Returns: [x, y, z, vx, vy, vz]
```

## 2. ICRS (RA, Dec) ↔ Galactic

```python
from_icrs = agama.fromICRStoGalactic  # precomputed rotation matrix

# Positions only
l, b = agama.transformCelestialCoords(from_icrs, ra*d2r, dec*d2r)

# With proper motions
l, b, pml, pmb = agama.transformCelestialCoords(from_icrs,
    ra*d2r, dec*d2r, pmra, pmdec)

# With uncertainties (covariance propagation)
l, b, pml, pmb, sig_l, sig_b, corr = agama.transformCelestialCoords(from_icrs,
    ra*d2r, dec*d2r, pmra, pmdec,
    sig_ra, sig_dec, corr_radec)
```

## 3. Cartesian ↔ Celestial (lon, lat, distance)

```python
# 3D → spherical
lon, lat, dist = agama.getCelestialCoords(x, y, z)
lon, lat, dist, pmlon, pmlat, vlos = agama.getCelestialCoords(x, y, z, vx, vy, vz)

# Spherical → 3D
x, y, z = agama.getCartesianCoords(lon, lat, dist)
x, y, z, vx, vy, vz = agama.getCartesianCoords(lon, lat, dist, pmlon, pmlat, vlos)
```

## 4. Rotation Matrices & Euler Angles

```python
# Create rotation matrix (intrinsic → observer)
R = agama.makeRotationMatrix(alpha, beta, gamma)  # Euler angles

# Extract Euler angles from matrix
alpha, beta, gamma = agama.getEulerAngles(R)

# Custom celestial pole rotation
R = agama.makeCelestialRotationMatrix(lon_pole, lat_pole, psi)
```

## 5. Practical Example: Star Coordinates

```python
# Convert a star's observed data to Galactocentric coordinates
# Input: RA=83.633°, Dec=+22.015°, plx=1.0 mas → d=1 kpc
# pmRA=20 mas/yr, pmDec=-5 mas/yr, vlos=10 km/s

ra, dec = 83.633, 22.015  # degrees
d = 1.0 / (1.0 * np.pi/648000)  # parallax → kpc (simplified)
pmra, pmdec = 20.0, -5.0  # mas/yr
vlos = 10.0  # km/s

# RA,Dec → Galactic
from_icrs = agama.fromICRStoGalactic
d2r = np.pi/180
l, b, d_kpc, pml, pmb, vlos_g = agama.transformCelestialCoords(
    from_icrs, ra*d2r, dec*d2r, pmra, pmdec)

# Galactic → Galactocentric
pos = agama.getGalactocentricFromGalactic(
    l, b, d_kpc, pml, pmb, vlos_g)
x, y, z, vx, vy, vz = pos