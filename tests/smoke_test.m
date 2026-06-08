% smoke_test.m
%
% Lightweight end-to-end check for the phase-drift workflow.
% This is not a validation sweep; it just verifies that the core routines
% run together on a very small problem and return finite diagnostics.

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

params