# Tests

This folder contains lightweight checks for the MATLAB workflow.

## Smoke test

Run from the repository root:

```matlab
run('tests/smoke_test.m')
```

The smoke test is intentionally small. It builds a short traveling-wave case, runs the phase-drift routine, and checks that the main diagnostics are present and finite. It is meant to catch basic path or interface problems, not to replace the validation sweeps in `experiments/`.
