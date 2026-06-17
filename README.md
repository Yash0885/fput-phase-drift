![phase drift signal](figures/phase_shift_vs_time.png)

Example phase shift signal extracted from the simulation.  
The approximately linear growth indicates a small propagation-speed mismatch.

# Phase drift in long-time symplectic simulation of an FPUT-beta traveling wave

This repository contains a small reproducible numerical experiment studying phase drift in long-time symplectic simulations of a traveling wave in the periodic Fermi-Pasta-Ulam-Tsingou (FPUT-beta) lattice.

The main question is why the numerical error appears to grow over long time intervals. In this experiment, two explanations are separated:

1. genuine deformation of the waveform  
2. accumulated phase drift caused by a small propagation-speed mismatch

The numerical evidence suggests that most of the observed error growth in this setup is caused by phase drift rather than visible deformation of the waveform.

---

# Results at a glance

The main diagnostic is to compare the numerical displacement profile against the reference wave before and after optimal spatial alignment.

| Diagnostic | Observed behavior |
|---|---|
| Direct waveform error | Grows over long time because the wave drifts relative to the fixed reference |
| Alignment-based waveform error | Remains much smaller after correcting for translation |
| Phase shift | Grows approximately linearly in time |
| Drift estimate | Small, typically on the order of `1e-5` in the tested runs |
| Energy drift | Remains very small under Störmer-Verlet time integration |

The interpretation is that the dominant long-time error is a coherent phase error, not an immediate loss of waveform shape.

---

# Model

We study the periodic FPUT-beta lattice with bond potential

```text
phi(r) = 0.5 r^2 + 0.25 r^4
```

and force

```text
phi'(r) = r + r^3
```

The equations of motion are

```text
x_ddot_i = phi'(x_{i+1} - x_i) - phi'(x_i - x_{i-1})
```

with periodic boundary conditions.

---

# Traveling wave construction

A traveling-wave-like profile is computed using

- Fourier spectral discretization
- Newton continuation
- a discrete traveling-wave residual equation

Parameters used in the main experiment:

```text
L = 16
k = pi / 16
speed shift target = 0.02
```

The reference wave speed is

```text
c_ref = c0 + speed_shift
```

where

```text
c0 = sqrt(2*(1 - cos(k))) / k
```

The base wave profile is repeated several times to create the full lattice used for time integration.

---

# Time integration

The system is integrated using the velocity-Verlet, or Störmer-Verlet, symplectic method.

Parameters used across the repo include:

```text
smoke-test dt = 0.01
time-step validation dt = 5e-3, 2.5e-3, 1.25e-3, 6.25e-4
simulation length about 1000 wave periods for validation runs
sampling rate = 4 samples per period
```

---

# Diagnostics

Several diagnostics are recorded during the simulation.

**Direct waveform error**  
Relative L2 difference between the numerical displacement profile and the fixed reference traveling wave profile.

**Alignment-based error**  
For each sampled displacement profile, the optimal spatial shift `s(t)` is computed by minimizing

```text
|| u(x,t) - u0(x + s) ||_2
```

Sub-grid translations are implemented using FFT interpolation.

**Phase drift estimate**  
The shift signal is

1. unwrapped,
2. corrected by subtracting the expected translation,
3. fit using linear regression.

The fitted slope gives a drift-rate estimate in this alignment-shift convention. Its sign depends on the convention used for the spatial shift, so it should be interpreted as the measured coherent phase-drift rate rather than a separately derived wave speed.

**Energy drift**  
Energy is computed as

```text
H = sum(0.5 * v_i^2) + sum(phi(x_{i+1} - x_i))
```

and remains nearly conserved throughout the simulation.

---

# Repository structure

```text
src/
    build_traveling_wave_frac.m
    run_phase_drift_from_tw_frac.m
    nsoli.m

experiments/
    validate_time_step.m
    validate_repeats.m
    validate_newton_tolerance.m
    generate_figures.m

tests/
    smoke_test.m

figures/
    generated plots

archive/
    earlier exploratory scripts
```

The `src` directory contains reusable functions.  
The `experiments` directory contains scripts that run the numerical tests.  
The `tests` directory contains lightweight checks for path and interface problems.  
The `figures` directory stores generated plots.

---

# Reproducing the experiment

From the repository root, run the following commands in MATLAB:

```matlab
addpath('src')
run('experiments/validate_time_step.m')
run('experiments/validate_repeats.m')
run('experiments/validate_newton_tolerance.m')
run('experiments/generate_figures.m')
```

The validation scripts write `.mat` and `.csv` result files. The figure-generation script reuses saved results if they already exist.

The expected output figures are:

```text
figures/direct_vs_aligned_error.png
figures/phase_shift_vs_time.png
figures/energy_drift.png
figures/drift_vs_dt.png
figures/drift_vs_repeats.png
figures/drift_vs_newton_tol.png
```

---

# Smoke test

A small smoke test is included for checking that the main MATLAB functions still run together on a short problem:

```matlab
run('tests/smoke_test.m')
```

This test is meant to catch basic path or interface issues. It is not a replacement for the longer validation sweeps in `experiments`.

---

# Notes and limitations

This repository documents a reproducible numerical experiment rather than a polished software package. The emphasis is on clear numerical diagnostics and interpretation of the results.

The results should be read as evidence for this controlled FPUT-beta traveling-wave setup, not as a general theorem for all nonlinear lattices, amplitudes, or numerical methods. The main point is that translation alignment separates phase drift from waveform deformation in a useful and measurable way.
