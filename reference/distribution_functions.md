# Distribution Functions Reference

## Overview
DistributionFunctions describe the probability distribution of stars in phase space (position + velocity). They are combined with a Potential via GalaxyModel to sample particles or compute moments.

```python
df = agama.DistributionFunction(params_dict)
# Composite DF
df = agama.DistributionFunction(df1, df2)
```

## 1. QuasiSpherical (type='quasispherical')

Eddington inversion: f(E) for a spherical density in a spherical potential, optionally with anisotropy.

```python
# Simplest: isotropic, potential provides density
df = agama.DistributionFunction(type='quasispherical', potential=pot)

# With explicit density (for composite potentials where dens ≠ pot)
df = agama.DistributionFunction(type='quasispherical',
    potential=pot, density=dens)

# With anisotropy
df = agama.DistributionFunction(type='quasispherical',
    potential=pot, density=dens,
    beta=0.3,              # constant anisotropy (0=iso, <0=tangential, >0=radial)
    anisotropyRadius=10.0, # Osipkov-Merritt: beta(r)=r²/(r²+ra²)
    rotFrac=0.5,           # fraction of stars rotating (0-1)
    Jphi0=1000)            # characteristic angular momentum for rotation
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `potential` | — | Potential object (required) |
| `density` | =pot | Density object (can be same as potential) |
| `beta` | 0 | Constant anisotropy parameter |
| `anisotropyRadius` | ∞ | Osipkov-Merritt anisotropy radius |
| `rotFrac` | 0 | Rotation fraction (0=non-rotating, 1=max rotation) |
| `Jphi0` | 0 | Angular momentum scale for rotation |

## 2. DoublePowerLaw (type='doublepowerlaw')

Analytic action-based DF: f(J) = norm × [1 + (J₀/h)ᶰ]^(-p) × exp[-(h/J_c)ᵠ]

where h(J) = (coefJrJz·Jr + coefJzJphi·|Jz| + (J_phi - Jphi0)) and J_phi = Lz.

```python
df = agama.DistributionFunction(type='doublepowerlaw',
    mass=1e12, J0=1000, slopeIn=1.5, slopeOut=3.5, steepness=2.0)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mass` | — | Total mass (alternative to norm) |
| `norm` | — | Normalization |
| `J0` | 1 | Break action |
| `Jcutoff` | ∞ | Cutoff action |
| `slopeIn` / `p` | 0.5 | Inner slope |
| `slopeOut` | 5 | Outer slope |
| `steepness` / `η` | 2 | Transition steepness |
| `coefJrJz` | 1 | Jr/Jz weighting in h(J) |
| `coefJzJphi` / `coefJzLz` | 0 | |Jz|/Lz weighting |
| `Jphi0` | 0 | Reference Lz for rotation |
| `rotFrac` | 0 | Rotation fraction |
| `Jn` | 0 | Additional slope parameter |

## 3. Exponential (type='exponential')

Simple action-based: f(J) ∝ exp(-h(J)/J₀). Same parameter structure as DoublePowerLaw.

```python
df = agama.DistributionFunction(type='exponential',
    mass=1e10, J0=500, coefJrJz=1.5)
```

## 4. QuasiIsothermal (type='quasiisothermal')

For disk populations: f(Jr, Jz, Lz) = Σ/σr² × exp(-Jr/σr - Jz/σz)

```python
df = agama.DistributionFunction(type='quasiisothermal',
    Sigma0=1e9, sigma0=30, sigma1=15, Rsigma=3.0)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Sigma0` | — | Central surface density |
| `mass` | — | Total mass (alternative to Sigma0) |
| `sigma0` | 1 | Radial velocity dispersion at R=0 |
| `sigma1` | 0 | Slope of sigma(R): sigma(R) = sigma0 × exp(-R/sigma1) |
| `Rsigma` | 1 | Scale length for sigma(R) |
| `R_sigma` | ∞ | Cutoff radius for sigma(R) |

## 5. Key Differences

| DF Type | Best For | Requires ActionFinder? | Speed |
|---------|----------|----------------------|-------|
| `quasispherical` | Spheroids, halos, bulges | No (uses E) | Fast |
| `quasiisothermal` | Disks | No (uses approximate actions) | Medium |
| `doublepowerlaw` | General, action-based | Yes (needs action finder) | Slow |
| `exponential` | Simple action-based | Yes | Slow-medium |