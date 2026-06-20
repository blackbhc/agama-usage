# GalaxyModel & Orbit Integration Reference

## 1. GalaxyModel

Combines Potential + DistributionFunction to sample particles, compute moments, and project kinematics.

### Constructor
```python
gm = agama.GalaxyModel(potential, df)
# potential: agama.Potential object
# df: agama.DistributionFunction object
```

### Sampling N-body Initial Conditions
```python
xv, masses = gm.sample(N)
# xv:    (N, 6) array [x, y, z, vx, vy, vz]
# masses: (N,) array, sum ≈ total mass of DF

# For component-by-component sampling (composite galaxy):
# Each component has its own Potential + DF
disk_gm  = agama.GalaxyModel(pot_disk,  df_disk)
halo_gm  = agama.GalaxyModel(pot_halo,  df_halo)
bulge_gm = agama.GalaxyModel(pot_bulge, df_bulge)

xv_disk, m_disk   = disk_gm.sample(N_disk)
xv_halo, m_halo   = halo_gm.sample(N_halo)
xv_bulge, m_bulge = bulge_gm.sample(N_bulge)

# Concatenate
xv_all   = np.vstack([xv_disk, xv_halo, xv_bulge])
m_all    = np.concatenate([m_disk, m_halo, m_bulge])
```

### Computing Velocity Moments
```python
xyz = np.random.randn(1000, 3) * 5  # positions
mom = gm.moments(xyz,
    dens=True,   # density at each point
    vel=True,    # mean velocity (3 components)
    vel2=True)   # velocity dispersion tensor (6 components: σxx, σyy, σzz, σxy, σxz, σyz)
# mom: (N, K) array
# dens -> column 0
# vel  -> columns 1-3: [vx, vy, vz]
# vel2 -> columns 4-9: [σxx, σyy, σzz, σxy, σxz, σyz]

# Example: get just density + radial velocity dispersion
rho_vr2 = gm.moments(xyz, dens=True, vel2=True)
```

### Line-of-Sight Velocity Distribution
```python
xy = np.random.randn(100, 2) * 3        # projected positions
vlos = np.linspace(-300, 300, 101)      # km/s
vdf = gm.vdf(xy, vlos, los_direction=2) # LOS along z-axis (2 = z)
```

## 2. Orbit Integration: agama.orbit()

### Basic Single Orbit
```python
times, traj = agama.orbit(
    potential=pot,
    ic=[x, y, z, vx, vy, vz],      # initial condition
    time=10.0,                      # integration time (internal units)
    trajsize=1001)                  # number of output snapshots
# times: (1001,) array
# traj:  (1001, 6) array [x, y, z, vx, vy, vz]
```

### Orbit Object (adaptive internal storage)
```python
orbit = agama.orbit(
    potential=pot,
    ic=[x, y, z, vx, vy, vz],
    time=10.0,
    dtype=object)      # returns Orbit object instead of arrays

# Access trajectory
full_traj = orbit(orbit)              # at internal steps
interp = orbit(np.linspace(0, 10, 100))  # interpolate at specific times
endpoint = orbit.end                  # final state
```

### Rotating (Barred) Potentials
```python
# Pattern speed in km/s/kpc, negative = clockwise
Omega_bar = -40.0

times, traj = agama.orbit(
    potential=pot,
    ic=[x, y, z, vx, vy, vz],
    time=10.0,
    Omega=Omega_bar,          # rotating frame
    trajsize=10000)

# Positions are in the co-rotating frame
# For visualization, plot in bar frame:
x_bar = traj[:, 0]
y_bar = traj[:, 1]
```

### Batch Orbits (Multiple ICs)
```python
ic_array = np.random.randn(20, 6) * np.array([3, 3, 1, 50, 50, 30])
times_arr, traj_arr = agama.orbit(
    potential=pot,
    ic=ic_array,              # (N, 6) array
    time=10.0,
    trajsize=501)

# times_arr: (N, 501) or (501,) depending on separateTime param
# traj_arr:  (N, 501, 6)

# With separate output per particle:
# times_arr: (N, 501)
# traj_arr:  (N, 501, 6)
```

### Full Parameter List

| Parameter | Default | Description |
|-----------|---------|-------------|
| `potential` | — | Potential object (required) |
| `ic` | — | Initial condition: (6,) or (N,6) array |
| `time` | — | Integration time (required) |
| `trajsize` | 0 | Number of output snapshots (0=default adaptive) |
| `dtype` | float | Output type: float for arrays, `object` for Orbit |
| `Omega` | 0 | Pattern speed (rotating frame) |
| `timestart` | 0 | Start time |
| `separateTime` | True | Separate time arrays per particle |
| `der` | False | Also return derivative vectors |
| `lyapunov` | False | Compute Lyapunov exponent |
| `targets` | None | List of Target objects for orbit storage |
| `accuracy` | 1e-8 | Integration accuracy |
| `maxStep` | 0 | Max timestep |
| `minStep` | 0 | Min timestep |
| `maxInterval` | 0 | Max output interval |

### Time-Dependent Potentials
```python
# Evolving potential from a list of snapshots
evolving = agama.Potential(type='Evolving',
    potential=[pot1, pot2, pot3],  # list of Potential objects
    time=[0, 5, 10],              # snapshot times
    interpLinear=False)            # False=piecewise, True=linear

times, traj = agama.orbit(
    potential=evolving,
    ic=ic,
    time=10.0,
    trajsize=1000)
```

### Lyapunov Exponents & Chaos Detection
```python
time, endpoint, dev_vectors, (lyap, chaos_time) = agama.orbit(
    potential=pot,
    ic=ic,
    time=100.,
    trajsize=1,       # only keep endpoint
    der=True,         # return deviation vectors
    lyapunov=True)    # compute Lyapunov exponent
# lyap: estimated Lyapunov exponent
# chaos_time: time when orbit is classified as chaotic
```

### Orbit Properties
```python
# Circular velocity at radius R
vc = np.sqrt(-R * pot.force(R, 0, 0)[0])

# Circular period at radius R
Tcirc = pot.Tcirc(R)

# Total energy
E = pot.potential(x, y, z) + 0.5 * (vx**2 + vy**2 + vz**2)

# Plot orbit in 2D projection
import matplotlib.pyplot as plt
plt.plot(traj[:, 0], traj[:, 1], lw=0.3)
plt.xlabel('x [kpc]'); plt.ylabel('y [kpc]')
plt.axis('equal')
```

## 3. ActionFinder

```python
af = agama.ActionFinder(pot)

# Actions only
J = af(ic_array)           # (N, 3): [Jr, Jz, Jphi]

# With frequencies
J, freqs = af(ic_array, frequencies=True)   # [Ωr, Ωz, Ωphi]

# With angles
J, angles = af(ic_array, angles=True)       # [θr, θz, θphi]

# All together
J, freqs, angles = af(ic_array, frequencies=True, angles=True)
```

| Param | Default | Description |
|-------|---------|-------------|
| `interp` | False | Use interpolation (faster but approximate) |

## 4. Self-Consistent Modeling Workflow

```python
# 1. Initial potential guess
pot = agama.Potential(...)
# 2. DF for each component
df_disk = agama.DistributionFunction(type='quasiisothermal', ...)
df_halo = agama.DistributionFunction(type='quasispherical', ...)
# 3. Sample particles
gm_disk = agama.GalaxyModel(pot, df_disk)
xv_disk, m_disk = gm_disk.sample(N_disk)
# 4. Create new potential from particle distribution
pot_new = agama.Potential(type='Multipole', particles=xv_disk, ...)
# 5. Iterate until convergence (compare density/potential)