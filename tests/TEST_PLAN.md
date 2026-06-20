# Agama Skill Test Plan

Three test scenarios to validate the agama-usage skill.

## Directory Structure
```
tests/
├── TEST_PLAN.md              # This file
├── test1_exponential_disk_nfw/
│   ├── plan.md               # Step-by-step plan
│   ├── run_test1.py          # Main script
│   └── results/              # Output plots & data
├── test2_satellite_merger/
│   ├── plan.md
│   ├── run_test2.py
│   └── results/
└── test3_bar_orbits/
    ├── plan.md
    ├── run_test3.py
    └── results/
```

## Quick Verification
Before any test: `python -c "import agama; agama.setUnits(1,1,1); print('OK')"`

## Validation Criteria
- Each script runs without errors
- Output files are created (plots, data files)
- Physical quantities are reasonable (rotation curves, orbital periods)