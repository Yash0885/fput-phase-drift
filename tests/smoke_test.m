% smoke_test.m
%
% Lightweight end-to-end check for the phase-drift workflow.
% This is not a validation sweep; it only checks that the core routines
% run together on a small problem and return finite diagnostics.

clear;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo_root, 'src'));

params.sh_target = 1e-4;
params.L         = 8;
params.k         = pi/8;
params.n_periods = 1;
params.n_repeats = 2;
params.dx        = 1.0;
params.dt        = 1e-2;
params.M         = 2;
params.align_tol = 1e-8;

params.phip = @(s) s + s.^3;
params.Phi  = @(s) 0.5*s.^2 + 0.25*s.^4;

params.newton_tol   = [1e-10, 1e-8];
params.newton_parms = [20, 20, 0.9, 3, 10];

tw = build_traveling_wave_frac(params);
out = run_phase_drift_from_tw_frac(tw, params);

required_fields = {'delta_c_driftOnly', 'E_align_max', ...
    'E_direct_max', 'energy_rel_max', 't_hist', 'shift_dist_driftOnly'};

for j = 1:numel(required_fields)
    assert(isfield(out, required_fields{j}), ...
        'Missing output field: %s', required_fields{j});
end

assert(numel(out.t_hist) >= 2, 'Smoke test produced too few samples.');
assert(all(isfinite(out.t_hist)), 'Nonfinite time history values.');
assert(all(isfinite(out.shift_dist_driftOnly)), 'Nonfinite shift values.');
assert(isfinite(out.E_align_max), 'Nonfinite aligned error.');
assert(isfinite(out.E_direct_max), 'Nonfinite direct error.');
assert(isfinite(out.energy_rel_max), 'Nonfinite energy drift.');

fprintf('smoke_test passed: %d samples, max aligned error %.3e\n', ...
    numel(out.t_hist), out.E_align_max);
