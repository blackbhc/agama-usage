# agama-usage

> **AGAMA (Action-based GAlaxy Modeling Architecture)** — a comprehensive C++/Python library for galactic dynamics: potentials, distribution functions, N-body initial conditions, orbit integration, action/angle computation, and self-consistent modeling.

**Version:** 1.0  
**Author:** Based on AGAMA by Eugene Vasiliev  
**Triggers:** `agama`, `N-body`, `galactic dynamics`, `galaxy`, `potential`, `orbit integration`, `distribution function`, `action-angle`, `tidal stream`, `self-consistent model`

## Table of Contents

1. [Quick Start](#1-quick-start)
2. [Units](#2-units)
3. [Density & Potential Models](#3-density--potential-models)
4. [Distribution Functions](#4-distribution-functions)
5. [GalaxyModel: Sampling N-body ICs & Moments](#5-galaxymodel)
6. [Orbit Integration](#6-orbit-integration)
7. [Action/Angle Variables](#7-actionangle-variables)
8. [Time-Dependent Potentials](#8-time-dependent-potentials)
9. [Coordinate Systems](#9-coordinate-systems)
10. [N-body & Test Particle Simulations](#10-n-body--test-particle-simulations)
11. [I/O & Utilities](#11-io--utilities)
12. [Complete Worked Examples](#12-complete-worked-examples)
13. [Reference](#13-reference)

---

## 1. Quick Start

```python
import agama
import numpy as np

# Set units: length=kpc, velocity=km/s, mass=Msun
agama.setUnits(length=1, velocity=1, mass=1)
# Time unit = length/velocity ≈ 0.978 Gyr per (kpc/(km/s))
# i.e., 1 time unit = 1 kpc / (1 km/s) ≈ 0.978 Gyr
TIME_UNIT_GYR = 0.978  # multiply time by this to get Gyr

# Create a simple Milky-Way like galaxy
disk = {'type': 'Disk', 'mass': 5e10, 'scaleRadius': 3.0, 'scaleHeight': 0.3}
halo = {'type': 'Spheroid', 'mass': 1e12, 'scaleRadius': 15.0, 'gamma': 1, 'beta': 3}
pot = agama.Potential(disk, halo)

# Compute circular velocity
R = np.linspace(0.1, 30, 100)
xyz = np.column_stack([R, np.zeros_like(R), np.zeros_like(R)])
vc = np.sqrt(-R * pot.force(xyz)[:, 0])

# Create isotropic DF and sample particles
df = agama.DistributionFunction(type='quasispherical', potential=pot)
gm = agama.GalaxyModel(pot, df)
xv, masses = gm.sample(10000)

# Integrate an orbit
times, traj = agama.orbit(potential=pot, ic=xv[0], time=10.0, trajsize=1001)
```

## 2. Units

### 2.1 Setting Units

```python
agama.setUnits(length=1, velocity=1, mass=1)     # kpc, km/s, Msun (default)
agama.setUnits(length=8, velocity=220, mass=1e10)  # or whatever you prefer
```

**Time is derived**: `time_unit = length_unit / velocity_unit`.  
With defaults: `1 unit = 1 kpc / (1 km/s) ≈ 0.978 Gyr`.

The `agama.G` constant is set automatically from units:  
`G ≈ 4.302 × 10⁻³ (km/s)² kpc / Msun` in default units.

### 2.2 Astropy Support

```python
import astropy.units as u
agama.setUnits(length=8*u.kpc, velocity=220*u.km/u.s, mass=1e10*u.Msun)
```

### 2.3 Unit Conversion Functions

```python
agama.getUnits()  # returns dict with astropy Quantities if available
```

## 3. Density & Potential Models

All density/potential models are created through `agama.Density(params)` or `agama.Potential(params)` where `params` is a dict with `type` and type-specific parameters.

### 3.1 Common Parameters

| Parameter | Description | Applies To |
|-----------|-------------|-----------|
| `type` | Model name (see below) | All |
| `mass` | Total mass | Most analytic models |
| `scaleRadius` / `rscale` | Scale/break radius | Most |
| `file` | Path to INI file with saved coefficients | All |
| `center` | Shift center: `[x,y,z]` | Modifier |
| `orientation` | Euler angles: `[α,β,γ]` | Modifier |
| `rotation` | Time-dependent rotation: `[[t0,ϕ0],[t1,ϕ1]]` | Modifier |
| `scale` | Rescale density/potential by factor | Modifier |

### 3.2 Density-Only Models

#### Disk (type='Disk')
Formula: `ρ(R,z) = f(R) × h(z)` where  
`f(R) = Σ₀ exp[-(R/R_d)^(1/n) - R_0/R + ε cos(R/R_d)]`  
`h(z) = exp(-|z/h|)` for h>0, `sech²(z/(2h))` for h<0, `δ(z)` for h=0

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `surfaceDensity` / `Sigma0` | — | mass/length² | Surface density at R=0 (alternative to `mass`) |
| `mass` | — | mass | Total mass (alternative to `surfaceDensity`) |
| `scaleRadius` / `rscale` | 1 | length | Disk scale length |
| `scaleHeight` / `scaleRadius2` | 1 | length | Scale height (negative → sech² profile) |
| `innerCutoffRadius` | 0 | length | Inner hole radius |
| `sersicIndex` | 1 | — | Sérsic index (n=1 = exponential) |
| `modulationAmplitude` | 0 | — | Amplitude of cos(R/R_d) modulation |

#### Spheroid (type='Spheroid')
Formula: `ρ(r) = ρ₀ (r/r₀)^(-γ) [1 + (r/r₀)^α]^((γ-β)/α) exp[-(r/r_cut)^ξ]`  
(NB: axis ratios p=y/x, q=z/x applied as ellipsoidal radius)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `densityNorm` / `rho0` | — | Density normalization |
| `mass` | — | Total mass (alternative to `densityNorm`) |
| `scaleRadius` / `rscale` | 1 | Break radius r₀ |
| `gamma` | 1 | Inner power-law slope (γ < 3) |
| `beta` | 4 | Outer power-law slope |
| `alpha` | 1 | Sharpness of transition |
| `outerCutoffRadius` | ∞ | Exponential cutoff radius |
| `cutoffStrength` / `xi` | 2 | Cutoff steepness |
| `axisRatioY` / `p` | 1 | Y/X axis ratio |
| `axisRatioZ` / `q` | 1 | Z/X axis ratio |

Special cases:
- **NFW**: `gamma=1, beta=3` (but use `type='NFW'` for exact NFW potential)
- **Einasto**: `alpha=1/n, gamma=0, beta=∞, cutoffStrength=1/n`
- **Dehnen**: `alpha=1, beta=4, gamma=[0,1,2,3]` (but use `type='Dehnen'` for exact)

#### Nuker (type='Nuker')
Deprojected Nuker surface brightness profile.

| Parameter | Description |
|-----------|-------------|
| `surfaceDensity` / `Sigma0` | Surface density at scale radius |
| `mass` | Total mass (alternative) |
| `scaleRadius` | Break radius |
| `gamma` | Inner slope (0≤γ<2) |
| `beta` | Outer slope (β>2) |
| `alpha` | Transition sharpness |
| `outerCutoffRadius` | Cutoff radius |
| `axisRatioY`, `axisRatioZ` | Axis ratios |

#### Sersic (type='Sersic')
Deprojected Sérsic profile: `Σ(R) = Σ₀ exp[-b (R/R_e)^(1/n)]`.

| Parameter | Description |
|-----------|-------------|
| `surfaceDensity` / `Sigma0` | Central surface density |
| `mass` | Total mass (alternative) |
| `scaleRadius` | Effective radius R_e |
| `sersicIndex` / `n` | Sérsic index |
| `axisRatioY`, `axisRatioZ` | Axis ratios |

### 3.3 Analytic Potential Models (also provide density)

#### Plummer (type='Plummer')
`Φ(r) = -GM / √(r² + b²)`

| Parameter | Description |
|-----------|-------------|
| `mass` | Total mass |
| `scaleRadius` / `rscale` | Scale radius b |

#### NFW (type='NFW')
Navarro-Frenk-White: `ρ(r) = ρ₀ / [(r/a)(1 + r/a)²]`

| Parameter | Description |
|-----------|-------------|
| `mass` | Alternative to `densityNorm` |
| `densityNorm` / `rho0` | Characteristic density |
| `scaleRadius` / `rscale` | Scale radius a |

The NFW potential is: `Φ(r) = -4πGρ₀ a² ln(1 + r/a) / (r/a)`

#### Dehnen (type='Dehnen')
`ρ(r) = (3-γ)M / (4πa³) (r/a)^(-γ) (1 + r/a)^(-(4-γ))`

| Parameter | Description |
|-----------|-------------|
| `mass` | Total mass |
| `scaleRadius` / `rscale` | Scale radius a |
| `gamma` | Inner slope (0 ≤ γ ≤ 3; γ=1 = Hernquist, γ=2 = Jaffe) |
| `axisRatioY`, `axisRatioZ` | Axis ratios |

#### Ferrers (type='Ferrers')
`ρ ∝ [1 - (r/a)²]^n` for r ≤ a, 0 outside. n=2 by default.

| Parameter | Description |
|-----------|-------------|
| `mass` | Total mass |
| `scaleRadius` / `rscale` | Radius a |
| `axisRatioY`, `axisRatioZ` | Axis ratios |

Note: Can be used as a simple bar model (triaxial Ferrers).

#### Miyamoto-Nagai (type='MiyamotoNagai')
`Φ(R,z) = -GM / √[R² + (a + √(z² + b²))²]`

| Parameter | Description |
|-----------|-------------|
| `mass` | Total mass |
| `scaleRadius` / `rscale` | Radial scale a |
| `scaleHeight` / `scaleRadius2` | Vertical scale b |

#### Isochrone (type='Isochrone')
`Φ(r) = -GM / [b + √(r² + b²)]`

| Parameter | Description |
|-----------|-------------|
| `mass` | Total mass |
| `scaleRadius` / `rscale` | Scale radius b |

#### King (type='King')
King model (lowered isothermal sphere).

| Parameter | Description |
|-----------|-------------|
| `mass` | Total mass |
| `scaleRadius` / `rscale` | Core radius |
| `W0` | Central potential depth |
| `trunc` | Truncation factor (default=1) |

#### Logarithmic (type='Logarithmic')
`Φ(r) = ½ v₀² ln[r_c² + x² + (y/p)² + (z/q)²]`

| Parameter | Description |
|-----------|-------------|
| `v0` | Asymptotic circular velocity |
| `scaleRadius` / `rscale` | Core radius r_c |
| `axisRatioY` / `p` | Y/X axis ratio |
| `axisRatioZ` / `q` | Z/X axis ratio |

#### Harmonic (type='Harmonic')
`Φ(r) = ½ Ω² [x² + (y/p)² + (z/q)²]`

| Parameter | Description |
|-----------|-------------|
| `Omega` | Oscillation frequency |
| `axisRatioY`, `axisRatioZ` | Axis ratios |

#### Perfect Ellipsoid (type='PerfectEllipsoid')
Homeoidal density distribution.

| Parameter | Description |
|-----------|-------------|
| `mass` | Total mass |
| `scaleRadius` | Scale radius |
| `axisRatioZ` / `q` | Z/X axis ratio |

#### KeplerBinary (type='KeplerBinary')
Time-dependent binary black hole.

| Parameter | Description |
|-----------|-------------|
| `mass` | Total binary mass |
| `binary_q` / `q` | Mass ratio (0 ≤ q ≤ 1) |
| `binary_sma` / `sma` | Semimajor axis |
| `binary_ecc` / `ecc` | Eccentricity |
| `binary_phase` / `phase` | Initial orbital phase |

### 3.4 Potential Expansions

#### Multipole (type='Multipole')
Spherical-harmonic expansion. Good for spheroidal/flattened systems.

```python
# From analytic density
pot = agama.Potential(type='Multipole', density=my_density_model,
                       gridsizer=30, rmin=0.01, rmax=100, lmax=6, mmax=6)
# From N-body particles
pot = agama.Potential(type='Multipole', particles=posvel_array,
                       gridsizer=30, rmin=0.01, rmax=100, lmax=8)
# From another potential
pot = agama.Potential(type='Multipole', potential=pot_original,
                       gridsizer=30, rmin=0.01, rmax=100, lmax=6)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `density` / `potential` / `particles` | — | Source: Density, Potential, or (N,6) array |
| `gridSizeR` | 25 | Radial grid nodes |
| `rmin` | auto | Minimum radius |
| `rmax` | auto | Maximum radius |
| `lmax` | 6 | Angular order |
| `mmax` | 6 | Azimuthal order |
| `symmetry` | auto | Symmetry type (S/A/T/B/R/N) |
| `smoothing` | 0 | Smoothing parameter |

#### CylSpline (type='CylSpline')
Azimuthal Fourier expansion + 2D B-spline in (R,z). Best for disks and bars.

```python
pot = agama.Potential(type='CylSpline', density=disk_density,
                       mmax=8, gridsizeR=25, gridsizez=25,
                       Rmin=0.1, Rmax=40, zmin=0.05, zmax=20)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `density` / `potential` / `particles` | — | Source |
| `gridSizeR` | 25 | Radial grid |
| `gridSizez` | 25 | Vertical grid |
| `Rmin` | auto | Minimum R |
| `Rmax` | auto | Maximum R |
| `zmin` | auto | Minimum z |
| `zmax` | auto | Maximum z |
| `mmax` | 6 | Azimuthal Fourier order |
| `symmetry` | auto | Symmetry type |

#### BasisSet (type='BasisSet')
SCF-style radial + spherical harmonic expansion.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nmax` | 6 | Radial order |
| `lmax` | 6 | Angular order |
| `mmax` | 6 | Azimuthal order |
| `r0` | 1 | Radial scale |
| `eta` | 0 | Radial basis parameter |

### 3.5 Creating Composite Potentials

```python
# From multiple dicts
pot = agama.Potential(
    {'type': 'Spheroid', 'mass': 2e10, 'scaleRadius': 0.5, 'gamma': 0, 'beta': 1.8},
    {'type': 'Disk',     'mass': 5e10, 'scaleRadius': 3.0, 'scaleHeight': -0.4},
    {'type': 'Spheroid', 'mass': 1e12, 'scaleRadius': 15.0, 'gamma': 1, 'beta': 3}
)

# From multiple Potential objects
pot1 = agama.Potential(type='Plummer', mass=1e6, scaleRadius=0.01)
pot2 = agama.Potential(type='NFW', mass=1e12, scaleRadius=15)
composite = agama.Potential(pot1, pot2)

# Index components
component0 = pot[0]  # first component
```

### 3.6 Reading/Writing Potentials to File

```python
# Save
pot = agama.Potential(type='Multipole', density=my_den, ...)
pot.export('my_potential.ini')  # or agama.writePotential('my_potential.ini', pot)

# Load
pot = agama.Potential(file='my_potential.ini')
```

## 4. Distribution Functions

### 4.1 QuasiSpherical (type='quasispherical')
Eddington inversion: f(E) for a given spherical density+potential, optionally with anisotropy.

```python
# Most common: isotropic DF from potential+density
df = agama.DistributionFunction(type='quasispherical', potential=pot, density=dens)
# If density is omitted, potential is used for both
df = agama.DistributionFunction(type='quasispherical', potential=pot)

# With anisotropy
df = agama.DistributionFunction(type='quasispherical',
    potential=pot, density=dens,
    beta=0.3,              # constant anisotropy parameter
    anisotropyRadius=10.0, # Osipkov-Merritt anisotropy radius
    rotFrac=0.5,           # rotation fraction
    Jphi0=1000)            # angular momentum for rotation
```

### 4.2 DoublePowerLaw (type='doublepowerlaw')
Analytic f(J) = norm × [1 + (J₀/h(J))^η]^(-p) × exp[-(h(J)/J_cut)^q]

| Parameter | Default | Description |
|-----------|---------|-------------|
| `norm` | — | Normalization |
| `mass` | — | Total mass (alternative to norm) |
| `J0` | 1 | Break action |
| `Jcutoff` | ∞ | Cutoff action |
| `slopeIn` | 0.5 | Inner slope (p) |
| `slopeOut` | 5 | Outer slope |
| `steepness` | 2 | Transition steepness (η) |
| `coefJrJz` | 1 | Ratio Jᵣ/J_z in h(J) |
| `coefJzJphi` / `coefJzLz` | 0 | |z|/J_φ in h(J) |
| `Jphi0` | 0 | Reference Lz |
| `rotFrac` | 0 | Rotation fraction |

### 4.3 Exponential (type='exponential')
f(J) ∝ exp(-h(J)/J₀). Same parameter structure as DoublePowerLaw.

### 4.4 QuasiIsothermal (type='quasiisothermal')
For disk populations: f(Jᵣ, J_z, Lz) = Σ/σᵣ² × exp(-Jᵣ/σᵣ - J_z/σ_z)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Sigma0` | — | Surface density normalization |
| `mass` | — | Total mass (alternative) |
| `sigma0` | 1 | Velocity dispersion scale |
| `sigma1` | 0 | Slope of σ(R) |
| `Rsigma` | 1 | Scale for σ(R) |
| `R_sigma` | ∞ | Cutoff for σ(R) |

### 4.5 Composite DF

```python
df = agama.DistributionFunction(df1, df2)
```

## 5. GalaxyModel

`GalaxyModel` combines a `Potential` and a `DistributionFunction` to compute velocity moments, sample particles, or compute projected kinematics.

### 5.1 Sampling N-body Initial Conditions

```python
gm = agama.GalaxyModel(potential, df)
xv, masses = gm.sample(N)  # N-body ICs
# xv: (N, 6) array [x, y, z, vx, vy, vz]
# masses: (N,) array, normalized so sum ≈ total mass
```

The sampling uses adaptive rejection sampling. For large N, it can be slow; for quick tests ~10⁴ particles it works well.

### 5.2 Velocity Moments

```python
xyz = np.random.randn(100, 3) * 5  # 100 points
mom = gm.moments(xyz, dens=True, vel=True, vel2=True)
# Returns: (N, K) array where K depends on requested outputs
# dens=True  → column 0: density
# vel=True   → columns 1-3: mean velocity (vx, vy, vz)
# vel2=True  → columns 4-9: velocity dispersion tensor (symmetric 3x3)
# Or use gm.moments(xyz, dens=False, vel=False, vel2=True) for just velocity dispersion
```

### 5.3 Projected DF (LOSVD)

```python
xy = np.random.randn(100, 2) * 3  # projected positions
vlos_grid = np.linspace(-300, 300, 101)  # km/s
vdf = gm.vdf(xy, vlos_grid, los_direction=2)  # LOS along z
```

### 5.4 Target Function (for Schwarzschild modeling)

```python
# Create a density target (surface density on a grid)
target = agama.Target(type='density', grid_parameters...)
stored = gm(target)  # store orbit contributions
```

## 6. Orbit Integration

### 6.1 Basic Usage

```python
# Single orbit — return timestamps + trajectory
times, traj = agama.orbit(potential=pot,
                           ic=[x, y, z, vx, vy, vz],
                           time=10.0,   # integration time
                           trajsize=1001)  # number of output points
# times: (1001,) array
# traj:  (1001, 6) array [x, y, z, vx, vy, vz]

# As an Orbit object (adaptive timestep, stores full internal state)
orbit = agama.orbit(potential=pot, ic=ic, time=10.0, dtype=object)
traj = orbit(orbit)           # get stored trajectory at full resolution
traj_interp = orbit(times)    # interpolate at specific timestamps
```

### 6.2 Rotating (Barred) Potentials

```python
# Omega sets the pattern speed (clockwise if negative)
times, traj_bar = agama.orbit(potential=pot, ic=ic, time=10.0,
                               Omega=40.0,   # km/s/kpc
                               trajsize=1001)
# Positions are in the rotating frame
```

### 6.3 Batch Orbits

```python
# Multiple orbits at once (much faster than looping)
ic_array = np.random.randn(100, 6) * np.array([3, 3, 1, 50, 50, 30])
times_arr, traj_arr = agama.orbit(potential=pot, ic=ic_array,
                                   time=10.0, trajsize=501)
# times_arr: (100, 501) or (501,) depending on separateTime
# traj_arr:  (100, 501, 6)
```

### 6.4 Orbit Properties

```python
# Circular velocity at radius
vc = np.sqrt(-R * pot_mw.force(R, 0, 0)[0])

# Circular orbit actions/energy (rough estimate)
Tcirc = pot.Tcirc(R)  # circular period at radius R

# Total energy
E = pot.potential(x, y, z) + 0.5 * (vx**2 + vy**2 + vz**2)
```

### 6.5 Time-Dependent Potentials

```python
# Create evolving potential
evol_pot = agama.Potential(type='Evolving',
    file='sequence_of_potentials.ini')  # or pass dict/list directly
times, traj = agama.orbit(potential=evol_pot, ic=ic, time=10.0, trajsize=1001)
```

## 7. Action/Angle Variables

### 7.1 ActionFinder

```python
# Create action finder for a potential
af = agama.ActionFinder(pot)

# Compute actions at positions
J = af(ic_array)  # (N, 3) array: [Jr, Jz, Jphi]

# With frequencies
J, freqs = af(ic_array, frequencies=True)  # [Ωr, Ωz, Ωphi]

# With angles
J, angles = af(ic_array, angles=True)

# All together: actions, frequencies, angles
J, fr, ang = af(ic_array, frequencies=True, angles=True)
```

### 7.2 Action Finder Types

| Parameter | Description |
|-----------|-------------|
| `interp=False` | Exact Staeckel Fudge (slower, more accurate) |
| `interp=True` | Interpolated (faster, approximate) |

## 8. Time-Dependent Potentials

### 8.1 Evolving Potential

```python
# From list of snapshots
evolving = agama.Potential(type='Evolving',
    potential=[pot1, pot2, pot3],
    time=[0, 5, 10],  # snapshot times
    interpLinear=False)  # False = piecewise constant, True = linear
```

### 8.2 Modifiers

Any potential can be modified:

```python
# Moving center
pot_moving = agama.Potential(type='Plummer',
    mass=1e6, scaleRadius=0.5,
    center=[[0, 0, 0, 0], [5, 5, 0, 10]])  # [x, y, z] at times [0, 10]

# Rotating (time-dependent orientation)
pot_rotating = agama.Potential(type='Ferrers',
    mass=1e10, scaleRadius=3, axisRatioY=0.3,
    rotation=[[0, 0], [10, 360*np.pi/180]])  # [time, angle]

# Scaling
pot_scaled = agama.Potential(type='Plummer',
    mass=1e6, scaleRadius=0.5,
    scale=[[0, 1], [5, 0.5], [10, 0.2]])  # [time, scale_factor]
```

## 9. Coordinate Systems

### 9.1 Galactocentric ↔ Galactic (l, b)

```python
# Galactocentric → Galactic
l, b, d, pml, pmb, vlos = agama.getGalacticFromGalactocentric(
    x, y, z, vx, vy, vz,
    galcen_distance=8.122,
    galcen_v_sun=(12.9, 245.6, 7.78),  # U, V, W
    z_sun=0.0208)

# Galactic → Galactocentric
pos = agama.getGalactocentricFromGalactic(
    l, b, d, pml, pmb, vlos,
    galcen_distance=8.122,
    galcen_v_sun=(12.9, 245.6, 7.78),
    z_sun=0.0208)
```

### 9.2 ICRS (RA, Dec) ↔ Galactic

```python
# Precomputed rotation matrix
from_icrs = agama.fromICRStoGalactic

# Transform
l, b = agama.transformCelestialCoords(from_icrs, ra * np.pi/180, dec * np.pi/180)

# With proper motions
l, b, pml, pmb = agama.transformCelestialCoords(from_icrs,
    ra*np.pi/180, dec*np.pi/180, pmra, pmdec)

# With covariance
l, b, pml, pmb, sig_l, sig_b, corr = agama.transformCelestialCoords(from_icrs,
    ra*d2r, dec*d2r, pmra, pmdec, sig_ra, sig_dec, corr_radec)
```

### 9.3 Cartesian ↔ Celestial

```python
lon, lat, dist = agama.getCelestialCoords(x, y, z)
lon, lat, dist, pmlon, pmlat, vlos = agama.getCelestialCoords(x, y, z, vx, vy, vz)
x, y, z = agama.getCartesianCoords(lon, lat, dist)
x, y, z, vx, vy, vz = agama.getCartesianCoords(lon, lat, dist, pmlon, pmlat, vlos)
```

### 9.4 Euler Angles & Rotation

```python
R = agama.makeRotationMatrix(alpha, beta, gamma)  # intrinsic → observer
alpha, beta, gamma = agama.getEulerAngles(R)

# For celestial coordinates with user-defined pole
R = agama.makeCelestialRotationMatrix(lon_pole, lat_pole, psi)
```

## 10. N-body & Test Particle Simulations

### 10.1 Simple N-body Initial Conditions

```python
# 1. Create a composite galaxy model
disk_pot = {'type': 'Disk', 'mass': 5e10, 'scaleRadius': 3, 'scaleHeight': 0.3}
halo_pot = {'type': 'Spheroid', 'mass': 1e12, 'scaleRadius': 15, 'gamma': 1, 'beta': 3}
pot = agama.Potential(disk_pot, halo_pot)
df = agama.DistributionFunction(type='quasispherical', potential=pot)
gm = agama.GalaxyModel(pot, df)
xv, masses = gm.sample(50000)

# 2. Write to file
agama.writeSnapshot('ICs.gadget', (xv, masses), 'gadget')
```

### 10.2 Satellite + Host Galaxy (Test Particle Simulation)

This simulates a satellite galaxy orbiting in a host potential. The satellite is represented by test particles that feel both the host potential and the satellite's own potential.

```python
# Host potential
pot_host = agama.Potential(
    {'type': 'Spheroid', 'mass': 2e10, 'scaleRadius': 0.5, 'gamma': 0, 'beta': 1.8},
    {'type': 'Disk',     'mass': 5e10, 'scaleRadius': 3, 'scaleHeight': -0.4},
    {'type': 'Spheroid', 'mass': 1e12, 'scaleRadius': 15, 'gamma': 1, 'beta': 3})

# Satellite potential + DF
pot_sat_init = agama.Potential(type='Plummer', mass=1e9, scaleRadius=0.5)
df_sat = agama.DistributionFunction(type='quasispherical', potential=pot_sat_init)
gm_sat = agama.GalaxyModel(pot_sat_init, df_sat)
xv_sat, m_sat = gm_sat.sample(1000)

# Place satellite at initial orbit position
R0 = 50.0  # kpc
V0 = 0.7 * np.sqrt(-R0 * pot_host.force(R0, 0, 0)[0])  # 70% of circular
xv_sat += np.array([R0, 0, 0, 0, V0, 0])

# Integrate particles as they orbit in host potential (test particle approx)
# For each particle, the force comes only from host + satellite's own potential
# placed at the satellite center (which moves on its own predefined orbit)
```

For a more detailed restricted N-body simulation (including dynamical friction, mass loss, and tidal stream), see the [Complete Worked Examples](#124-satellite-tidal-stream) section.

### 10.3 Dynamical Friction (Chandrasekhar)

```python
def dynamical_friction(pos, vel, satellite_mass, host_pot, sigma_fn):
    """Compute Chandrasekhar dynamical friction acceleration."""
    dist = np.sqrt(np.sum(pos**2))
    v_mag = np.sqrt(np.sum(vel**2))
    rho = host_pot.density(pos)
    sigma = sigma_fn(dist)
    X = v_mag / (sigma * np.sqrt(2))
    import scipy.special
    coulomb_log = np.log(dist / 2.0)  # approximate
    return -(4 * np.pi * rho * vel / v_mag *
        (scipy.special.erf(X) - 2/np.pi**0.5 * X * np.exp(-X**2)) *
        satellite_mass * agama.G**2 / v_mag**2 * coulomb_log)
```

## 11. I/O & Utilities

### 11.1 Snapshot Formats

```python
# Read Gadget format
pos, vel, mass, header = agama.readSnapshot('snapshot_001.gadget')
# pos: (N, 3), vel: (N, 3), mass: (N,)

# Write Gadget format
agama.writeSnapshot('out.gadget', (xv, masses), 'gadget')
# or using combined pos+vel array
agama.writeSnapshot('out.gadget', (pos_vel_array, masses), 'gadget')
```

### 11.2 Spline Utilities (from pygama.py)

```python
# Create non-uniform grid
grid = agama.nonuniformGrid(30, 0.1, 100)  # 30 nodes from ~0.1 to 100
grid = agama.symmetricGrid(30, 0.1, 50)    # symmetric about 0

# Interpolation with AGAMA's cubic spline
spl = agama.Spline(x_grid, y_values)
y_interp = spl(x_new)
```

### 11.3 Utility Functions

```python
from agama import Potential, Density, DistributionFunction, GalaxyModel
from agama import readSnapshot, writeSnapshot, setUnits, getUnits
from agama import fromICRStoGalactic, transformCelestialCoords
from agama import getGalacticFromGalactocentric, getGalactocentricFromGalactic
from agama import getCelestialCoords, getCartesianCoords
from agama import nonuniformGrid, symmetricGrid
from agama import makeRotationMatrix, getEulerAngles
from agama.G = 4.302e-3  # already set by setUnits in (km/s)² kpc / Msun
```

## 12. Complete Worked Examples

### 12.1 Galaxy with Exponential Disk + NFW Halo

```python
import agama, numpy as np, matplotlib.pyplot as plt
agama.setUnits(length=1, velocity=1, mass=1)  # kpc, km/s, Msun

# Create galaxy
disk = {'type': 'Disk', 'mass': 5e10, 'scaleRadius': 3.0, 'scaleHeight': 0.3}
halo = {'type': 'Spheroid', 'mass': 1e12, 'scaleRadius': 15.0, 'gamma': 1, 'beta': 3}
pot = agama.Potential(disk, halo)

# Circular velocity
R = np.linspace(0.1, 30, 100)
xyz = np.column_stack([R, np.zeros_like(R), np.zeros_like(R)])
vc = np.sqrt(-R * pot.force(xyz)[:, 0])
plt.plot(R, vc); plt.xlabel('R [kpc]'); plt.ylabel('Vc [km/s]')
plt.savefig('vcirc.pdf')

# Sample N-body ICs
df = agama.DistributionFunction(type='quasispherical', potential=pot)
gm = agama.GalaxyModel(pot, df)
xv, masses = gm.sample(10000)

# Plot particle distribution
plt.figure()
plt.scatter(xv[:, 0], xv[:, 1], s=0.5, alpha=0.5)
plt.xlabel('x [kpc]'); plt.ylabel('y [kpc]')
plt.savefig('particles.pdf')
```

### 12.2 Orbits in a Barred Potential

```python
import agama, numpy as np, matplotlib.pyplot as plt
agama.setUnits(length=1, velocity=1, mass=1)

# Create a bar using Ferrers ellipsoid + disk + halo
bar = {'type': 'Ferrers', 'mass': 2e10, 'scaleRadius': 3.0,
       'axisRatioY': 0.3, 'axisRatioZ': 0.2}
disk = {'type': 'Disk', 'mass': 5e10, 'scaleRadius': 3.0, 'scaleHeight': 0.3}
halo = {'type': 'Spheroid', 'mass': 1e12, 'scaleRadius': 15.0, 'gamma': 1, 'beta': 3}
pot = agama.Potential(bar, disk, halo)

# Pattern speed (km/s/kpc, negative = clockwise)
Omega_bar = -40.0

# 20 test particles with different initial conditions
np.random.seed(42)
ics = []
for i in range(20):
    R = 2.0 + np.random.rand() * 6.0  # 2-8 kpc
    phi = np.random.rand() * 2 * np.pi
    # Circular velocity estimate
    vc = np.sqrt(-R * pot.force(R, 0, 0)[0])
    vr = (np.random.rand() - 0.5) * 40  # km/s radial perturbation
    vt = vc + (np.random.rand() - 0.5) * 30  # km/s tangential perturbation
    x = R * np.cos(phi)
    y = R * np.sin(phi)
    vx = vr * np.cos(phi) - vt * np.sin(phi)
    vy = vr * np.sin(phi) + vt * np.cos(phi)
    ics.append([x, y, 0, vx, vy, 0])

ics = np.array(ics)

# Integrate
times, trajs = agama.orbit(potential=pot, ic=ics, time=10.0,
                           Omega=Omega_bar, trajsize=10000)

# Plot in rotating frame (positions are already in rotating frame)
fig, axes = plt.subplots(5, 4, figsize=(16, 20))
for i in range(5):
    for j in range(4):
        k = i * 4 + j
        ax = axes[i][j]
        if k < 20:
            ax.plot(trajs[k, :, 0], trajs[k, :, 1], lw=0.5)
            ax.set_xlim(-10, 10)
            ax.set_ylim(-10, 10)
            ax.set_aspect('equal')
        ax.set_title(f'Particle {k+1}')

plt.tight_layout()
plt.savefig('bar_orbits.pdf')
```

### 12.3 Action-Based N-body IC

```python
import agama, numpy as np
agama.setUnits(length=1, velocity=1, mass=1)

# Create potential
pot = agama.Potential(type='Multipole',
    density={'type': 'Spheroid', 'mass': 1e12, 'scaleRadius': 15, 'gamma': 1, 'beta': 3},
    gridsizer=30, rmin=0.1, rmax=200, lmax=6)

# Create DoublePowerLaw DF fitted to reproduce the density
df = agama.DistributionFunction(type='doublepowerlaw',
    mass=1e12, J0=1000, slopeIn=1.5, slopeOut=3.5, steepness=2.0)

# Sample
gm = agama.GalaxyModel(pot, df)
xv, masses = gm.sample(20000)
```

### 12.4 Satellite Tidal Stream

```python
import agama, numpy as np, matplotlib.pyplot as plt
agama.setUnits(length=1, velocity=1, mass=1)

# Host galaxy
pot_host = agama.Potential(
    {'type': 'Spheroid', 'mass': 2e10, 'scaleRadius': 0.5, 'gamma': 0, 'beta': 1.8},
    {'type': 'Disk', 'mass': 5e10, 'scaleRadius': 3.0, 'scaleHeight': 0.3},
    {'type': 'Spheroid', 'mass': 1e12, 'scaleRadius': 15.0, 'gamma': 1, 'beta': 3})

# Satellite (Plummer)
pot_sat = agama.Potential(type='Plummer', mass=1e9, scaleRadius=0.5)
df_sat = agama.DistributionFunction(type='quasispherical', potential=pot_sat)
gm_sat = agama.GalaxyModel(pot_sat, df_sat)
xv_sat, m_sat = gm_sat.sample(1000)

# Initial orbit (apocenter at 50 kpc, eccentricity ~0.3)
R0 = 50.0
Vcirc = np.sqrt(-R0 * pot_host.force(R0, 0, 0)[0])
V0 = 0.7 * Vcirc
xv_sat += np.array([R0, 0, 0, 0, V0, 0])

# Integrate satellite center of mass only (test particles in static host)
times, traj_sat = agama.orbit(potential=pot_host,
    ic=[R0, 0, 0, 0, V0, 0], time=6.0, trajsize=1000)

# For real test particle integration, one would need to integrate each particle
# in the time-dependent potential of host + moving satellite
```

## 13. Reference

### 13.1 Parameter Table Summary

| Python Parameter | C++ Parameter | Notes |
|-----------------|--------------|-------|
| `surfaceDensity` | `Sigma0` | Alternative names via `popString(..., "surfaceDensity", "Sigma0")` |
| `scaleRadius` | `rscale` | |
| `scaleHeight` | `scaleRadius2` | Negative → sech², positive → exp, zero → thin |
| `densityNorm` | `rho0` | |
| `axisRatioY` | `p` | |
| `axisRatioZ` | `q` | |
| `cutoffStrength` | `xi` | |
| `outerCutoffRadius` | — | Infinity by default |

### 13.2 Key Source Files

| File | Content |
|------|---------|
| `src/potential_factory.cpp` | All potential/density creation logic, parameter parsing |
| `src/potential_base.h` | Base classes for Density and Potential |
| `src/potential_analytic.h` | Plummer, NFW, Dehnen, Ferrers, Isochrone, etc. |
| `src/potential_spheroid.h` | Spheroid (alpha-beta-gamma), Nuker, Sersic |
| `src/potential_disk.h` | Exponential/disk density and potential |
| `src/df_factory.cpp` | DF creation (QuasiSpherical, DoublePowerLaw, etc.) |
| `src/galaxymodel_base.h` | GalaxyModel class |
| `src/orbit.h` | Orbit integration |
| `src/actions_staeckel.h` | Staeckel Fudge for actions |
| `src/interface_python.cpp` | All Python bindings |
| `py/pygama.py` | Python extras (coords, grids, wrappers) |

### 13.3 Units System

Default: `length=1 kpc, velocity=1 km/s, mass=1 Msun`

| Quantity | Expression | Value in default units |
|----------|-----------|----------------------|
| Time | `length / velocity` | 0.978 Gyr |
| G | gravitational constant | 4.302 × 10⁻³ |
| Action | `length × velocity` | 0.978 kpc·km/s·Gyr ≈ 1 kpc·km/s |

### 13.4 Error Handling

- Unknown parameters → `RuntimeError` with list of unknown keys
- Invalid parameter values → `RuntimeError` with description
- DF without potential → `invalid_argument`
- Expansions from particles → `RuntimeError` if too few particles

### 13.5 Performance Notes

- `Multipole` potential evaluation: O(lmax²) — use lmax ≤ 6-10
- `CylSpline` potential evaluation: O(mmax × gridResolution) — use mmax ≤ 12
- `ActionFinder` with `interp=False`: ~10x slower than `interp=True`
- Batch orbit integration is much faster than looping (uses vectorized internal operations)
- Sampling from GalaxyModel: adaptive rejection sampling can be slow for large N; use ~10⁴ for tests
- Default time unit: 1.0 in code ≈ 0.978 Gyr

### 13.6 References

- Vasiliev, E. (2019). AGAMA: Action-based galaxy modelling architecture. MNRAS, 482(2), 1525-1544.
- Dehnen, W. & Binney, J. (1998). Mass models of the Milky Way. MNRAS, 294(3), 429-438.
- Binney, J. & Tremaine, S. (2008). Galactic Dynamics (2nd ed.). Princeton University Press.
