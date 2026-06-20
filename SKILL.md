---
name: agama-usage
description: |
  AGAMA (Action-based GAlaxy Modeling Architecture) galactic dynamics library usage — potentials, distribution functions, N-body initial conditions, orbit integration, test particle simulations, action/angles, and self-consistent modeling.

  Use when: agama, N-body, galactic dynamics, galaxy modeling, potential construction, orbit integration, distribution function, tidal stream, action-angle, test particle simulation, self-consistent model
allowed-tools:
  - Read
  - Write
  - Bash
  - Edit
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

## Workflow

### Phase 1: Setup & Units
1. Call `agama.setUnits(length=1, velocity=1, mass=1)` at the start (default: kpc, km/s, Msun)
2. Remember: time unit = length/velocity (~0.978 Gyr for defaults)
3. Create potentials as dicts with `type`, `mass`, `scaleRadius`, etc.

### Phase 2: Build Galaxy Potential
Use `agama.Potential(...)` to create composite potentials from multiple components:
```python
pot = agama.Potential(
    {'type': 'Disk', 'mass': 5e10, 'scaleRadius': 3.0, 'scaleHeight': 0.3},
    {'type': 'NFW', 'mass': 1e12, 'scaleRadius': 15.0},
    {'type': 'Spheroid', 'mass': 2e10, 'scaleRadius': 0.5, 'gamma': 0, 'beta': 1.8})
```
See `reference/potentials.md` for all ~20 potential types and their parameters.

### Phase 3: Add Distribution Function & Sample
- For equilibrium ICs: create a `DistributionFunction` then a `GalaxyModel`
```python
df = agama.DistributionFunction(type='quasispherical', potential=pot)
gm = agama.GalaxyModel(pot, df)
xv, masses = gm.sample(N)  # (N,6) phase-space, (N,) masses
```
- For disk DF: use `type='quasiisothermal'` with `Sigma0`, `sigma0`, `sigma1`
- For action-based DF: use `type='doublepowerlaw'` or `type='exponential'`
See `reference/distribution_functions.md` for all DF types.

### Phase 4: Integrate Orbits
```python
times, traj = agama.orbit(potential=pot, ic=[x,y,z,vx,vy,vz],
                          time=10.0, trajsize=1001)
```
- For rotating (barred) potentials: pass `Omega=<pattern_speed>`
- For batch orbits: pass `ic` as (N,6) array
- For time-dependent potentials: use `agama.Potential(type='Evolving', ...)`

See `reference/galaxy_model.md` for detailed API.

### Phase 5: Actions & Angles
```python
af = agama.ActionFinder(pot)
J = af(ic)                              # (N,3) [Jr, Jz, Jphi]
J, freqs = af(ic, frequencies=True)     # + frequencies
J, angs = af(ic, angles=True)           # + angles
```

### Phase 6: Test Particle Simulation
1. Create host potential
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

## Examples Directory
Full worked examples are in `reference/examples.md`. Test scripts are in `examples/`.

## Complete Worked Examples (see reference/examples.md)
1. Exponential Disk + NFW Halo galaxy model
2. 1000-particle satellite merger test particles
3. Barred potential orbit integration (20 particles, 10 Gyr)

## Checklist
- [ ] Units set via `agama.setUnits()` before any other calls
- [ ] Time unit verified: `1 unit = length/velocity ≈ 0.978 Gyr` for defaults
- [ ] Potential parameters match documented types (check `reference/potentials.md`)
- [ ] QuasiSpherical DF provided with both `density` and `potential` if needed
- [ ] GalaxyModel constructed with correct Potential + DF combination
- [ ] Orbit `time` parameter in internal time units (not Gyr directly)
- [ ] For rotating potentials: `Omega` sign convention checked (negative = clockwise)
- [ ] Batch orbit integration uses (N,6) array, not list
- [ ] Snapshot I/O uses correct format string ('gadget', 'nemo')