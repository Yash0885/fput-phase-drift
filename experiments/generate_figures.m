addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'))

% generate_figures.m
%
% Generate the main figures for the phase-drift study.
%
% If a required results file is missing, run the corresponding experiment.

clear;
clc;

figdir = fullfile(pwd, 'figures');

if ~exist(figdir, 'dir')
    [ok, msg, msgID] = mkdir(figdir);
    if ~ok
        error('Could not create figures folder: %s (%s)', msg, msgID);
    end
end

%% ensure results exist

if ~isfile('validate_time_step_results.mat')
    fprintf('validate_time_step_results.mat not found -- running validate_time_step\n');
    validate_time_step
else
    fprintf('Loading existing validate_time_step_results.mat\n');
end

if ~isfile('validate_repeats_results.mat')
    fprintf('validate_repeats_results.mat not found -- running validate_repeats\n');
    validate_repeats
else
    fprintf('Loading existing validate_repeats_results.mat\n');
end

if ~isfile('validate_newton_tolerance_results.mat')
    fprintf('validate_newton_tolerance_results.mat not found -- running validate_newton_tolerance\n');
    validate_newton_tolerance
else
    fprintf('Loading existing validate_newton_tolerance_results.mat\n');
end

%% load results

S_dt  = load('validate_time_step_results.mat');
S_rep = load('validate_repeats_results.mat');
S_new = load('validate_newton_tolerance_results.mat');

outs_dt  = S_dt.outs;
outs_rep = S_rep.outs;
outs_new = S_new.outs;

%% representative run for time-history plots

[~, idx_fine] = min([outs_dt.dt]);
out = outs_dt(idx_fine);

fprintf('\nUsing dt = %.6g for time-history plots.\n', out.dt);

idx = 1:5:numel(out.t_hist);

%% direct vs aligned error

figure;

semilogy(out.t_hist(idx), out.E_direct(idx),'LineWidth',0.4); hold on;
semilogy(out.t_hist(idx), out.E_align(idx),'LineWidth',1.2);

grid on;

xlabel('time');
ylabel('relative error');

legend({'direct','aligned'},'Location','best');

title('Direct and aligned error vs time');

saveas(gcf, fullfile(figdir,'direct_vs_aligned_error.png'));

%% phase drift signal

figure;

plot(out.t_hist(idx), out.shift_dist_driftOnly(idx),'LineWidth',1.2); hold on;

j0 = max(2, round(0.1*numel(out.t_hist)));

p = polyfit(out.t_hist(j0:end), out.shift_dist_driftOnly(j0:end),1);

plot(out.t_hist, polyval(p,out.t_hist),'--','LineWidth',1.2);

grid on;

xlabel('time');
ylabel('shift distance');

legend({'drift signal','linear fit'},'Location','best');

title('Phase drift signal');

saveas(gcf, fullfile(figdir,'phase_shift_vs_time.png'));

%% energy drift

figure;

semilogy(out.t_hist(idx), out.Erel_hist(idx),'LineWidth',1.0);

grid on;

xlabel('time');
ylabel('relative energy drift');

title('Energy drift vs time');

saveas(gcf, fullfile(figdir,'energy_drift.png'));

%% drift estimate vs time step

figure;

dt_vals = [outs_dt.dt];
dc_vals = [outs_dt.delta_c_driftOnly];
dc_se   = [outs_dt.delta_c_se];

errorbar(dt_vals, dc_vals, dc_se,'o-','LineWidth',1.2);

grid on;

xlabel('time step dt');
ylabel('drift estimate');

title('Drift estimate vs time step');

saveas(gcf, fullfile(figdir,'drift_vs_dt.png'));

%% drift estimate vs repeated domain size

figure;

rep_vals = [outs_rep.n_repeats];
dc_rep   = 1e6 * [outs_rep.delta_c_driftOnly];
se_rep   = 1e6 * [outs_rep.delta_c_se];

errorbar(rep_vals, dc_rep, se_rep,'o-','LineWidth',1.2);

grid on;

xlabel('number of repeats');
ylabel('drift estimate x 10^6');

title('Drift estimate vs repeated-domain size');

saveas(gcf, fullfile(figdir,'drift_vs_repeats.png'));

%% drift estimate vs Newton tolerance

figure;

rtol_vals = [outs_new.rtol];
dc_new    = abs([outs_new.delta_c_driftOnly]);
se_new    = [outs_new.delta_c_se];

loglog(rtol_vals, dc_new,'o-','LineWidth',1.2); hold on;
loglog(rtol_vals, se_new,'s--','LineWidth',1.2);

grid on;

xlabel('relative Newton tolerance');
ylabel('magnitude');

legend({'absolute drift','standard error'},'Location','best');

title('Drift estimate vs Newton tolerance');

saveas(gcf, fullfile(figdir,'drift_vs_newton_tol.png'));

%% summary

fprintf('\nWrote figures to %s\n',figdir);
fprintf('  direct_vs_aligned_error.png\n');
fprintf('  phase_shift_vs_time.png\n');
fprintf('  energy_drift.png\n');
fprintf('  drift_vs_dt.png\n');
fprintf('  drift_vs_repeats.png\n');
fprintf('  drift_vs_newton_tol.png\n\n');
