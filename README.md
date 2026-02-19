# fput-phase-drift

Numerical experiments for phase drift in nonlinear traveling waves of the Fermi–Pasta–Ulam–Tsingou (FPUT) lattice.

The main diagnostic is:

* **best alignment shift vs time** (integer lattice sites), to measure slow phase drift
* **aligned shape error vs time**, to separate drift from true shape instability

## What’s in here

* `PhaseDrift.m` — single-case run that samples once per period and plots:

  * unwrapped best shift vs time
  * aligned relative error vs time
  * optional linear fit for drift rate

* `FPUT_repeat_sweep.m` — sweeps over `sh_target` values:

  * repeats the wave profile across a longer domain
  * integrates long time with Störmer–Verlet
  * reports final aligned error + phase shift
  * optional traces (shift/error/energy) for selected `sh_target`

## How to run

Open MATLAB, change directory into the repository folder, then run:

```matlab
addpath(genpath(pwd))
PhaseDrift
```

To run the parameter sweep:

```matlab
addpath(genpath(pwd))
FPUT_repeat_sweep
```

## Requirements

* MATLAB (tested on R2022+)
* No external toolboxes required

## Expected output

`PhaseDrift.m` should produce:

* Unwrapped best-shift (phase drift) vs time
* Aligned relative error vs time (log scale)

`FPUT_repeat_sweep.m` should produce:

* Stability plots over `sh_target` values
* Optional trace plots (shift/error/energy) for selected cases
