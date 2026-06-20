---
name: agama-usage
description: |
  AGAMA (Action-based GAlaxy Modeling Architecture) galactic dynamics library usage — potentials, distribution functions, N-body initial conditions, orbit integration, test particle simulations, action/angles, and self-consistent modeling.

  Use when: agama, N-body, galactic dynamics, galaxy modeling, potential construction, orbit integration, distribution function, tidal stream, action-angle, test particle simulation, self-consistent model
user-invocable: true
allowed-tools:
  - Read
  - Write
  - Bash
  - Edit
metadata:
  version: 1.1.0
  category: scientific-computing
---

# Skill: agama-usage

## Purpose
Guide an agent in using the AGAMA library for galactic dynamics: constructing galaxy models (potentials + distribution functions), generating N-body initial conditions, running test-particle simulations, integrating orbits (including in time-dependent/rotating potentials), computing action/angle variables, and performing self-consistent modeling.

## Use When
- User mentions "agama", "N-body", "galactic dynamics", "galaxy model", "potential", "orbit integration", "distribution function", "tidal stream", "action-angle", "self-consistent model"
- Task involves constructing a galaxy potential, sampling stellar particles, simulating satellite mergers, or integrating orbits in barred/spiral/evolving potentials

## Inputs
- Unit system (length, velocity, mass) — kpc, km/s, Msun by default
- Galaxy component parameters (mass, scale radius, profile shapes, axis ratios)
- Distribution function parameters (anisotropy, rotation fraction)
- Initial conditions for orbits or N-body particles
- Integration time and snapshot frequency

## Repository Structure
```
agama-usage/
  SKILL.md                  ← this file
  configs/                  ← .ini configuration files for galaxy models
    test1_exponential_disk_nfw.ini
    test2_mw_host.ini
    test3_bar_potential.ini
    McMillan17.ini           ← canonical MW model (from AGAMA data/)
    MWPotential2014.ini      ← Galpy MWPotential2014
    PriceWhelan17.ini        ← Gala MilkyWayPotential
    Cautun20.ini             ← Cautun+2020 MW with contracted halo
  reference/                ← detailed reference docs
    potentials.md
    distribution_functions.md
    galaxy_model.md
    coordinates.md
  examples/                 ← (planned) worked example scripts
```

## Workflow

### Phase 1: Setup & Units
1. Call `agama.setUnits(length=1, velocity=1, mass=1)` at the start (default: kpc, km/s, Msun)
2. Remember: time unit = length/velocity (~0.978 Gyr for defaults)
3. Set units BEFORE constructing any potential or DF objects

### Phase 2: Configuration via INI Files (Recommended)

AGAMA reads multi-component potentials from INI-format files. **Always prefer separate .ini files** for model parameters rather than hard-coding them in Python scripts. This keeps science configuration decoupled from analysis code.

#### INI File Format

Each section name must start with `Potential` (case-insensitive). Parameters are `key = value` pairs inside each section:

```ini
# comments start with # or ;
# Units: 1 kpc, 1 km/s, 1 Msun

[Potential disk]
type = Disk
mass = 5e10
scaleRadius = 3.0
scaleHeight = 0.3

[Potential halo]
type = NFW
mass = 1.2e12
scaleRadius = 15.0
```

Key rules:
- Section name MUST start with `Potential` (it can be `[Potential0]`, `[Potential disk]`, `[Potential bulge]`, etc.)
- Values are parsed as floats; scientific notation works (`5e10`)
- Comments with `#` at line start or inline
- Use `type` to specify the potential/density class
- Units must match what `agama.setUnits()` expects (default: kpc, km/s, Msun)

#### Loading from INI File

```python
# Single INI file with multiple [Potential*] sections → composite potential
pot = agama.Potential(file='path/to/model.ini')

# INI files can also be loaded via readPotential
# All sections starting with "Potential" are combined into one composite potential
```

#### When to Use INI vs Python Dict

| Approach | Best For |
|----------|----------|
| **INI file** (recommended) | Complex multi-component models; saving/reusing configurations; collaboration |
| **Python dict** | Quick tests; programmatic parameter generation; time-dependent modifiers |

#### Reference INI Files Included

This skill ships with pre-built reference INI files in `configs/`:

| File | Description |
|------|-------------|
| `McMillan17.ini` | McMillan 2017 best-fit MW model (thin+thick disk, gas, bulge, halo) |
| `MWPotential2014.ini` | Galpy MWPotential2014 (bulge + Miyamoto-Nagai disk + NFW halo) |
| `PriceWhelan17.ini` | Gala MilkyWayPotential (Miyamoto-Nagai disk + Dehnen bulge + NFW halo) |
| `Cautun20.ini` | Cautun+2020 MW with adiabatically contracted DM halo (uses Multipole) |
| `test1_exponential_disk_nfw.ini` | Simple exponential disk + NFW halo (for quick tests) |
| `test2_mw_host.ini` | Dehnen bulge + exponential disk + Logarithmic halo (flat rotation curve) |
| `test3_bar_potential.ini` | Bulge + disk + Ferrers bar + NFW halo |

These are also usable as starting points for custom models.

#### Advanced: INI with Time-Dependent Modifiers

```ini
[Potential bar]
type = Ferrers
mass = 2e10
scaleRadius = 4.0
axisRatioY = 0.35
# Time-dependent rotation (for bar pattern speed)
rotation = [[0, 0], [1, 0.7]]    # [time, angle] pairs (internal time units)

# Time-dependent center motion
center = [[0, 0, 0, 0], [10, 5, 0, 0]]    # [t, x, y, z]
```

#### INI with Expansion Coefficients

For complex potentials saved from a Multipole/CylSpline expansion:

```ini
[Potential halo]
type = Multipole
gridSizeR = 25
lmax = 0
symmetry = Spherical
# Coefficients block: tab-separated columns
# Phi
# radius    l=0,m=0
0.010000   -162153.09363897
0.017783   -162087.61429742
...
# dPhi/dr
# radius    l=0,m=0
0.010000   7593.92277406
0.017783   9108.19858139
...
```

Such files are produced by `pot.export('filename.ini')`.

### Phase 2 (alt): Build Galaxy Potential via Python Dict

When you need programmatic control (e.g., scanning parameters, time-dependent modifiers):

```python
pot = agama.Potential(
    {'type': 'Disk', 'mass': 5e10, 'scaleRadius': 3.0, 'scaleHeight': 0.3},
    {'type': 'NFW', 'mass': 1e12, 'scaleRadius': 15.0},
    {'type': 'Spheroid', 'mass': 2e10, 'scaleRadius': 0.5, 'gamma': 0, 'beta': 1.8})
```

### Phase 3: Add Distribution Function & Sample
- For equilibrium ICs: create a `DistributionFunction` then a `GalaxyModel`
```python
df = agama.DistributionFunction(type='quasispherical', potential=pot)
gm = agama.GalaxyModel(pot, df)
xv, masses = gm.sample(N)  # (N,6) phase-space, (N,) masses
```
- For disk DF: use `type='quasiisothermal'` with `Sigma0`, `sigma0`, `sigma1`
- For action-based DF: use `type='doublepowerlaw'` or `type='exponential'`

### Phase 4: Integrate Orbits
```python
times, traj = agama.orbit(potential=pot, ic=[x,y,z,vx,vy,vz],
                          time=10.0, trajsize=1001)
```
- For rotating (barred) potentials: pass `Omega=<pattern_speed>`
- For batch orbits: pass `ic` as (N,6) array

### Phase 5: Actions & Angles
```python
af = agama.ActionFinder(pot)
J = af(ic)                              # (N,3) [Jr, Jz, Jphi]
J, freqs = af(ic, frequencies=True)     # + frequencies
J, angs = af(ic, angles=True)           # + angles
```

### Phase 6: Test Particle Simulation
1. Create host potential (prefer INI file)
2. Create satellite potential + DF, sample particles
3. Offset satellite to initial orbital position
4. Integrate particles in time-dependent host+satellite potential
5. Track tidal stripping and orbital evolution

### Phase 7: I/O & Visualization
```python
agama.writeSnapshot('file.gadget', (xv, masses), 'gadget')
pos, vel, mass, hdr = agama.readSnapshot('file.gadget')
```

## Outputs
- Phase-space coordinates of N-body particles [(N,6) array]
- Orbit trajectories [(T,6) array per particle]
- Velocity moments at specified positions
- Distribution function values
- Action/angle variables
- Snapshot files in Gadget/NEMO format
- Plots of density profiles, rotation curves, orbital trajectories

## Limitations
- **Does NOT cover**: SPH/hydrodynamical simulations, radiative transfer, MHD, cosmological simulations
- **Does NOT provide**: GPU acceleration for AGAMA (CPU only)
- **Does NOT replace**: full N-body codes like GADGET/GyrfalcON for self-consistent N-body (AGAMA handles ICs + restricted N-body)
- **GalaxyModel.sample()** can be slow for N > 10^5; use ~10^4 for tests
- Time-dependent potentials with modifiers (`center`, `rotation`, `scale`) should be carefully tested

## Quick Reference

### Unit Conversions
| Task | Formula | Value (default units) |
|------|---------|-----------------------|
| Time | length/velocity | 1 unit ≈ 0.978 Gyr |
| G | constant | 4.302e-3 (km/s)² kpc/Msun |
| Action | length × velocity | 1 unit ≈ 1 kpc·km/s |

### Potential Types at a Glance
| type | Category | Key Params |
|------|----------|------------|
| `Disk` | Density | mass, scaleRadius, scaleHeight |
| `Spheroid` | Density | mass, scaleRadius, gamma, beta |
| `NFW` | Potential | mass/densityNorm, scaleRadius |
| `Plummer` | Potential | mass, scaleRadius |
| `Dehnen` | Potential | mass, scaleRadius, gamma |
| `Ferrers` | Potential | mass, scaleRadius, axisRatioY/Z |
| `MiyamotoNagai` | Potential | mass, scaleRadius, scaleHeight |
| `Logarithmic` | Potential | v0, scaleRadius |
| `Multipole` | Expansion | density/particles, lmax, gridsizeR |
| `CylSpline` | Expansion | density/particles, mmax, gridsizeR/Z |
| `BasisSet` | Expansion | density/particles, nmax, lmax |

### DF Types at a Glance
| type | When to Use | Key Params |
|------|-------------|------------|
| `quasispherical` | Spheroidal components | potential, density, beta, rotFrac |
| `quasiisothermal` | Disk components | Sigma0, sigma0, sigma1 |
| `doublepowerlaw` | Action-based general | J0, slopeIn, slopeOut |
| `exponential` | Simple action-based | J0 |

## Reference Docs
See `reference/` for detailed API references:
- `potentials.md` — all ~20 density/potential types with parameters
- `distribution_functions.md` — DF types, parameters, and usage
- `galaxy_model.md` — GalaxyModel, orbit integration, ActionFinder, SCM

## Examples Directory
Full worked examples are in `examples/` (planned). Reference INI configs are in `configs/`.

## Checklist
- [ ] Units set via `agama.setUnits()` before any other calls
- [ ] Time unit verified: `1 unit = length/velocity ≈ 0.978 Gyr` for defaults
- [ ] Potential parameters match documented types (check `reference/potentials.md`)
- [ ] For complex multi-component models, prefer separate `.ini` configuration files
- [ ] INI section names start with `Potential` (case-insensitive)
- [ ] QuasiSpherical DF provided with both `density` and `potential` if needed
- [ ] GalaxyModel constructed with correct Potential + DF combination
- [ ] Orbit `time` parameter in internal time units (not Gyr directly)
- [ ] For rotating potentials: `Omega` sign convention checked (negative = clockwise)
- [ ] Batch orbit integration uses (N,6) array, not list
- [ ] Snapshot I/O uses correct format string ('gadget', 'nemo')
