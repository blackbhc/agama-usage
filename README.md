# agama-usage

**AGAMA galactic dynamics skill — potentials, distribution functions, N-body ICs, orbit integration, test particle simulations, action/angles, and self-consistent modeling.**

This is a [Kun](https://kun.app) skill for guiding AI agents in using the [AGAMA](https://github.com/GalacticDynamics-Oxford/Agama) library (Action-based GAlaxy Modeling Architecture).

## Repository Structure

```
agama-usage/
├── SKILL.md                  # Main skill definition — workflow guide
├── README.md                 # This file
├── .gitignore
├── configs/                  # Pre-built .ini potential configuration files
│   ├── McMillan17.ini        # Canonical MW model (McMillan 2017)
│   ├── MWPotential2014.ini   # Galpy MWPotential2014
│   ├── PriceWhelan17.ini     # Gala MilkyWayPotential
│   ├── Cautun20.ini          # Cautun+2020 MW + contracted DM halo (Multipole)
│   ├── test1_exponential_disk_nfw.ini  # Simple disk+NFW template
│   ├── test2_mw_host.ini     # MW-like host (Dehnen+disk+Log halo)
│   └── test3_bar_potential.ini  # Bulge+disk+Ferrers bar+NFW halo
├── examples/                 # Standalone demo folders (each has script + config)
│   ├── demo1/           # Exp disk + NFW N-body ICs
│   ├── demo2/    # Binary BH test particle orbits
│   ├── demo3/    # Barred: box & loop orbits
│   └── demo4_elliptical_box_loop_orbits/  # Triaxial ellip: box & tube orbits
└── reference/                # Detailed API reference docs
    ├── potentials.md
    ├── distribution_functions.md
    ├── galaxy_model.md
    └── coordinates.md
```

## Key Principle

**Always prefer separate .ini configuration files** for model parameters rather than hard-coding them in Python scripts. This keeps science configuration decoupled from analysis code.

## Quick Start

```python
import agama

# 1. Set units (kpc, km/s, Msun)
agama.setUnits(length=1, velocity=1, mass=1)

# 2. Load a composite potential from INI
pot = agama.Potential(file='configs/test1_exponential_disk_nfw.ini')

# 3. Create a DF and sample particles
df = agama.DistributionFunction(type='quasispherical', potential=pot)
xv, masses = agama.GalaxyModel(pot, df).sample(10000)

# 4. Integrate an orbit
times, traj = agama.orbit(potential=pot, ic=[8, 0, 0, 0, 220, 0],
                          time=5.0, trajsize=1001)
```

## Contents

- **SKILL.md**: Complete workflow guide with INI configuration, potential types, DF types, orbit integration, and best practices.
- **configs/**: 7 reference `.ini` files for common MW models and test configurations.
- **reference/**: In-depth API docs for potentials, DFs, GalaxyModel, and coordinate systems.
- **examples/**: Runnable Python scripts demonstrating core workflows.

## Installation

AGAMA is a C++ library with Python bindings. Install from source:

```bash
git clone https://github.com/GalacticDynamics-Oxford/Agama.git
cd Agama
pip install .
```

Requires: numpy, scipy, matplotlib, GSL. For more details, refer to the official AGAMA repository: <https://github.com/GalacticDynamics-Oxford/Agama>

## License

MIT
