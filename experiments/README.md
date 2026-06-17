# Experiment scripts

This folder contains the longer validation runs used to reproduce the main phase-drift figures.

Run these scripts from the repository root in MATLAB:

```matlab
run('experiments/validate_time_step.m')
run('experiments/validate_repeats.m')
run('experiments/validate_newton_tolerance.m')
run('experiments/generate_figures.m')
```

## Scripts

- `validate_time_step.m` varies the time step while reusing the same traveling wave.
- `validate_repeats.m` varies the repeated-domain size while reusing the same base wave.
- `validate_newton_tolerance.m` rebuilds the traveling wave for several Newton tolerances.
- `generate_figures.m` loads the saved validation results, runs missing validation scripts if needed, and writes the main plots to `figures/`.

The validation scripts write `.mat` and `.csv` result files in the current working directory. For the intended layout, run them from the repository root.

These scripts are more expensive than the smoke test in `tests/`. Use `tests/smoke_test.m` first when checking basic path or interface problems.
