# Test 1: Exponential Disk + NFW Dark Matter Halo

**Reference**: arxiv:2602.19995

## Steps

### Step 1: Define Galaxy Parameters
- Units: kpc, km/s, Msun
- Disk: exponential, M_disk = 5×10¹⁰ Msun, R_d = 3.0 kpc, h_z = 0.3 kpc
- Halo: NFW, M_200 = 1.2×10¹² Msun, R_s = 15 kpc (concentration ~10)

### Step 2: Build Potentials
- Disk: `type='Disk'` with mass, scaleRadius, scaleHeight
- Halo: `type='NFW'` with mass (sets densityNorm internally), scaleRadius
- Combine into composite potential

### Step 3: Compute Rotation Curve
- V_c(R) = sqrt(-R * F_R) at R = [0.1, 30] kpc
- Plot: V_c vs R, decomposed by component (disk, halo, total)

### Step 4: Create Distribution Functions
- Halo: `type='quasispherical'` with the composite potential (isotropic)
- Disk: `type='quasiisothermal'` with Sigma0, sigma0=30, sigma1=15, Rsigma=3.0

### Step 5: Sample Particles
- Halo: 10000 particles via GalaxyModel
- Disk: 5000 particles via GalaxyModel
- Check: positions should be physically reasonable (R distribution, z distribution)

### Step 6: Density Profile Verification
- Bin particles radially, compare to input density profiles
- Plot: 3D density ρ(r) for halo, surface density Σ(R) for disk

### Step 7: Save ICs
- Save combined particle set to Gadget snapshot

### Step 8: Quick Orbit Check
- Integrate a test particle at R=8 kpc, circular velocity
- Verify it stays bound and orbits with roughly correct period

## Verification Criteria
- Rotation curve peaks at ~180-220 km/s
- Disk density profile follows Σ ∝ exp(-R/R_d)
- Halo density profile follows NFW shape
- Orbits are stable and bound