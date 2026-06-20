# Potential & Density Models Reference

## Construction Methods

### Method 1: INI Configuration Files (Recommended)
Keep model parameters in standalone `.ini` files:

```python
# Load all [Potential*] sections from an INI file as a composite potential
pot = agama.Potential(file='path/to/model.ini')
```

See **INI Format** below and the `configs/` directory in this skill repo for examples.

### Method 2: Python Dicts
```python
dens = agama.Density(params_dict)
pot  = agama.Potential(params_dict)
# Composite
pot = agama.Potential(dict1, dict2, ...)
pot = agama.Potential(pot_obj1, pot_obj2, ...)
```

## INI Format
```ini
# Lines starting with # are comments
# Section name must start with 'Potential' (case-insensitive)

[Potential anything here]
type = Spheroid         # required: density/potential class name
mass = 1e12             # model parameter
scaleRadius = 15.0      # another parameter
axisRatioZ = 0.8        # flattening

[Potential another component]
type = Disk
mass = 5e10
scaleRadius = 3.0
scaleHeight = 0.3
```

AGAMA ships with many INI files in its `data/` directory (McMillan17.ini, MWPotential2014.ini, etc.),
and this skill includes several reference INIs in `configs/`.

## 1. Common Modifier Parameters
## 1. Common Modifier Parameters

Applied to any potential as extra dict keys:

| Parameter | Description | Example |
|-----------|-------------|---------|
| `center` | Shift center [x,y,z] | `'center': [0,0,0]` |
| `orientation` | Euler angles [α,β,γ] | `'orientation': [0,0,0.3]` |
| `rotation` | Time-dep rotation [[t0,ϕ0],...] | `'rotation': [[0,0],[10,π]]` |
| `scale` | Time-dep scaling [[t,f],...] | `'scale': [[0,1],[5,0.5]]` |

## 2. Density-Only Models

### Disk (type='Disk')
```python
disk = {'type': 'Disk', 'mass': 5e10, 'scaleRadius': 3.0, 'scaleHeight': 0.3}
```
Formula: ρ(R,z) = Σ(R) × h(z)
- Σ(R) = Σ₀ exp[-(R/R_d)^(1/n) - R_0/R + ε cos(R/R_d)]
- h(z) = exp(-|z/h|) for h>0, sech²(z/(2h)) for h<0, δ(z) for h=0

| Parameter | Default | Description |
|-----------|---------|-------------|
| `surfaceDensity` / `Sigma0` | — | Surface density at R=0 |
| `mass` | — | Total mass (alternative) |
| `scaleRadius` / `rscale` | 1 | Scale length R_d |
| `scaleHeight` / `scaleRadius2` | 1 | Scale height (negative→sech², positive→exp, zero→thin) |
| `innerCutoffRadius` | 0 | Inner hole R_0 |
| `sersicIndex` | 1 | Sersic index (n=1 ≡ exponential) |
| `modulationAmplitude` | 0 | cos(R/R_d) modulation |

### Spheroid (type='Spheroid')
```python
spheroid = {'type': 'Spheroid', 'mass': 1e12, 'scaleRadius': 15,
            'gamma': 1, 'beta': 3, 'alpha': 1}
```
Formula: ρ(r) = ρ₀ (r/r₀)^(-γ) [1 + (r/r₀)^α]^((γ-β)/α) exp[-(r/r_cut)^ξ]

Ellipsoidal radius: r² = x² + (y/p)² + (z/q)²

| Parameter | Default | Description |
|-----------|---------|-------------|
| `densityNorm` / `rho0` | — | Density normalization |
| `mass` | — | Total mass (alternative) |
| `scaleRadius` / `rscale` | 1 | Break radius r₀ |
| `gamma` | 1 | Inner slope (γ < 3) |
| `beta` | 4 | Outer slope |
| `alpha` | 1 | Transition sharpness |
| `outerCutoffRadius` | ∞ | Exponential cutoff |
| `cutoffStrength` / `xi` | 2 | Cutoff steepness |
| `axisRatioY` / `p` | 1 | y/x ratio |
| `axisRatioZ` / `q` | 1 | z/x ratio |

Special cases mapped to Spheroid:
- **NFW-like**: gamma=1, beta=3 (but use `type='NFW'` for exact potential)
- **Einasto**: alpha=1/n, gamma=0, beta→∞, cutoffStrength=1/n
- **Plummer-like**: gamma=0, beta=5, alpha=2

### Nuker (type='Nuker')
Deprojected Nuker profile. Parameters: `surfaceDensity`/`mass`, `scaleRadius`, `gamma`, `beta`, `alpha`, `outerCutoffRadius`, `axisRatioY`, `axisRatioZ`.

### Sersic (type='Sersic')
Deprojected Sersic profile Σ(R) = Σ₀ exp[-b(R/R_e)^(1/n)].
Parameters: `surfaceDensity`/`mass`, `scaleRadius`(=R_e), `sersicIndex`/`n`, `axisRatioY`, `axisRatioZ`.

## 3. Analytic Potential Models

### Plummer (type='Plummer')
Φ(r) = -GM / √(r² + b²)

| Parameter | Description |
|-----------|-------------|
| `mass` | Total mass |
| `scaleRadius` / `rscale` | Plummer scale b |

### NFW (type='NFW')
ρ(r) = ρ₀ / [(r/a)(1 + r/a)²]
Φ(r) = -4πGρ₀ a² ln(1 + r/a) / (r/a)

| Parameter | Description |
|-----------|-------------|
| `mass` | Total mass (sets densityNorm internally) |
| `densityNorm` / `rho0` | Characteristic density ρ₀ |
| `scaleRadius` / `rscale` | Scale radius a |

### Dehnen (type='Dehnen')
ρ(r) = (3-γ)M/(4πa³) (r/a)^(-γ) (1+r/a)^(-(4-γ))

| Parameter | Description |
|-----------|-------------|
| `mass` | Total mass |
| `scaleRadius` / `rscale` | Scale radius a |
| `gamma` | Inner slope (0≤γ≤3; 1=Hernquist, 2=Jaffe) |
| `axisRatioY`, `axisRatioZ` | Axis ratios |

### Ferrers (type='Ferrers')
ρ ∝ [1 - (r/a)²]^n for r ≤ a, 0 outside. Default n=2.

| Parameter | Description |
|-----------|-------------|
| `mass` | Total mass |
| `scaleRadius` / `rscale` | Radius a |
| `axisRatioY`, `axisRatioZ` | Axis ratios (triaxial → simple bar) |

### Miyamoto-Nagai (type='MiyamotoNagai')
Φ(R,z) = -GM / √[R² + (a + √(z² + b²))²]

| Parameter | Description |
|-----------|-------------|
| `mass` | Total mass |
| `scaleRadius` / `rscale` | Radial scale a |
| `scaleHeight` / `scaleRadius2` | Vertical scale b |

### Long-Murali (type='LongMurali')
Barred generalization of Miyamoto-Nagai. Extra: `barLength`.

### Logarithmic (type='Logarithmic')
Φ = ½ v₀² ln[r_c² + x² + (y/p)² + (z/q)²]

| Parameter | Description |
|-----------|-------------|
| `v0` | Asymptotic circular velocity |
| `scaleRadius` / `rscale` | Core radius r_c |
| `axisRatioY` / `p`, `axisRatioZ` / `q` | Axis ratios |

Ideal for dark halo with flat rotation curve.

### Isochrone (type='Isochrone')
Φ(r) = -GM / [b + √(r² + b²)]

| Parameter | Description |
|-----------|-------------|
| `mass` | Total mass |
| `scaleRadius` / `rscale` | Scale b |

### King (type='King')
Lowered isothermal sphere with truncation.

| Parameter | Description |
|-----------|-------------|
| `mass` | Total mass |
| `scaleRadius` / `rscale` | Core radius |
| `W0` | Central potential depth |
| `trunc` | Truncation factor |

### Harmonic (type='Harmonic')
Φ = ½ Ω² [x² + (y/p)² + (z/q)²]

| Parameter | Description |
|-----------|-------------|
| `Omega` | Frequency |
| `axisRatioY`, `axisRatioZ` | Axis ratios |

## 4. Potential Expansions (from density/particles)

### Multipole (type='Multipole')
Spherical-harmonic expansion. Good for spheroidal/flattened systems.

```python
# From analytic density
pot = agama.Potential(type='Multipole',
    density=my_dens, gridsizeR=30, rmin=0.01, rmax=100, lmax=6)
# From particles
pot = agama.Potential(type='Multipole',
    particles=posvel_array, gridsizeR=30, rmin=0.01, rmax=100, lmax=8)
# From another potential
pot = agama.Potential(type='Multipole',
    potential=pot_original, gridsizeR=30, lmax=6)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `density`/`potential`/`particles` | — | Source |
| `gridSizeR` | 25 | Radial grid nodes |
| `rmin` | auto | Minimum radius |
| `rmax` | auto | Maximum radius |
| `lmax` | 6 | Angular order (cost O(l²)) |
| `mmax` | 6 | Azimuthal order |
| `symmetry` | auto | S/A/T/B/R/N |

### CylSpline (type='CylSpline')
Azimuthal Fourier + 2D B-spline in (R,z). Best for disks and bars.

```python
pot = agama.Potential(type='CylSpline',
    density=disk_dens, mmax=8, gridsizeR=25, gridsizez=25,
    Rmin=0.1, Rmax=40, zmin=0.05, zmax=20)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `gridSizeR` | 25 | Radial grid |
| `gridSizez` | 25 | Vertical grid |
| `Rmin`/`Rmax` | auto | Radial range |
| `zmin`/`zmax` | auto | Vertical range |
| `mmax` | 6 | Fourier order |

### BasisSet (type='BasisSet')
SCF-style radial + spherical harmonic expansion.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nmax` | 6 | Radial order |
| `lmax` | 6 | Angular order |
| `mmax` | 6 | Azimuthal order |
| `r0` | 1 | Radial scale |
| `eta` | 0 | Basis parameter |

## 5. Time-Dependent / Evolving Potentials

```python
# Evolving: interpolate between snapshots
evolving = agama.Potential(type='Evolving',
    potential=[pot1, pot2, pot3],
    time=[0, 5, 10],
    interpLinear=False)  # False=piecewise const, True=linear

# Modifiers for any potential
moving = agama.Potential(type='Plummer',
    mass=1e6, scaleRadius=0.5,
    center=[[0,0,0,0], [5,5,0,10]])  # [x,y,z] at times [0,10]
rotating = agama.Potential(type='Ferrers',
    mass=1e10, scaleRadius=3, axisRatioY=0.3,
    rotation=[[0,0],[10,360*np.pi/180]])  # [time, angle]
scaling = agama.Potential(type='Plummer',
    mass=1e6, scaleRadius=0.5,
    scale=[[0,1], [5,0.5], [10,0.2]])  # [time, factor]
```

## 6. Potential Queries

```python
pot = agama.Potential(...)
# Force (negative gradient)
Fx, Fy, Fz = pot.force(x, y, z)       # returns (3,) or (N,3)
# Density (if available)
rho = pot.density(x, y, z)
# Potential value
phi = pot.potential(x, y, z)
# Save/load
pot.export('pot.ini')
pot_loaded = agama.Potential(file='pot.ini')
# Components
comp = pot[0]  # first component of composite