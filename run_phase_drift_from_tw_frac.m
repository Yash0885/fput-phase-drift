function out = run_phase_drift_from_tw_frac(tw, params)
% run_phase_drift_from_tw_frac
%
% Evolve a supplied traveling wave and record drift/error diagnostics.

dt         = params.dt;
n_periods  = params.n_periods;
n_repeats  = params.n_repeats;
dx         = params.dx;
M          = params.M;
align_tol  = params.align_tol;
phip       = params.phip;
Phi        = params.Phi;

MAKE_PLOTS = false;
if isfield(params, 'MAKE_PLOTS')
    MAKE_PLOTS = params.MAKE_PLOTS;
end

uu    = tw.uu;
vper  = tw.vper;
c_ref = tw.c_ref;

N_wave = length(uu);
k = params.k;

% build the repeated domain from one base wave
x  = repmat(uu,   n_repeats, 1);
v  = repmat(vper, n_repeats, 1);
x0 = x;

% Fourier wave numbers for fractional shifts on one base period
q_wave = 2*pi * fftshift((-N_wave/2):(N_wave/2-1))' / N_wave;
q_wave = q_wave(:);

% effective stiffness proxy from the initial strains
strain0 = x0([2:end 1]) - x0;
kappa0  = phip_prime_vec_local(phip, strain0);

omega_eff0    = 2 * sqrt(max(kappa0));
omega_eff_max = omega_eff0;

% equilibrium linearization
phip_prime0 = estimate_phip_prime0_local(phip);
omega_max_eq = 2 * sqrt(max(phip_prime0, 0));

T_period = 2*pi / (c_ref * k);
steps    = round((n_periods * T_period) / dt);

% sample M times per wave period
store_every = max(1, round((T_period / dt) / M));
dt_sample   = store_every * dt;

nstore = floor(steps / store_every) + 1;

t_hist     = zeros(nstore, 1);
E_align    = zeros(nstore, 1);
E_direct   = zeros(nstore, 1);
shift_best = zeros(nstore, 1);
Erel_hist  = zeros(nstore, 1);

E0 = energy_FPU_local(x, v, Phi);
[~, a] = rhs_local(x, v, phip);

store_i = 1;
t_hist(store_i) = 0;
[E_align(store_i), shift_best(store_i)] = ...
    align_error_frac_local(x, x0, N_wave, q_wave, align_tol);
E_direct(store_i) = norm(x - x0) / max(norm(x0), 1);
Erel_hist(store_i) = 0;

% velocity-Verlet time stepping
for step = 1:steps

    x_new = x + v*dt + 0.5*dt^2*a;
    [~, a_new] = rhs_local(x_new, v, phip);
    v_new = v + 0.5*(a + a_new)*dt;

    x = x_new;
    v = v_new;
    a = a_new;

    if mod(step, store_every) == 0
        store_i = store_i + 1;

        t_hist(store_i) = step * dt;

        % best fractional shift against the base wave
        [E_align(store_i), shift_best(store_i)] = ...
            align_error_frac_local(x, x0, N_wave, q_wave, align_tol);

        E_direct(store_i) = norm(x - x0) / max(norm(x0), 1);

        E_now = energy_FPU_local(x, v, Phi);
        Erel_hist(store_i) = abs(E_now - E0) / max(1, abs(E0));

        strain = x([2:end 1]) - x;
        kappa  = phip_prime_vec_local(phip, strain);
        omega_eff_now = 2 * sqrt(max(kappa));
        omega_eff_max = max(omega_eff_max, omega_eff_now);
    end
end

t_hist     = t_hist(1:store_i);
E_align    = E_align(1:store_i);
E_direct   = E_direct(1:store_i);
shift_best = shift_best(1:store_i);
Erel_hist  = Erel_hist(1:store_i);

% build raw and driftOnly shift signals
if M <= 1
    shift_unwrap = unwrap_nearest_local(shift_best, N_wave);
    shift_dist_raw = shift_unwrap * dx;
    shift_dist_driftOnly = shift_dist_raw;
else
    % for M > 1 the sampled shift is only defined mod one base period
    % lock it near the expected translation before unwrapping
    m = (0:(numel(t_hist)-1))';
    expected_sites = -m * (c_ref * dt_sample / dx);
    expected_mod   = wrap_to_half_local(expected_sites, N_wave);

    phase_err_sites = wrap_to_half_local(shift_best - expected_mod, N_wave);

    phase_err_unwrap = unwrap_nearest_local(phase_err_sites, N_wave);
    shift_dist_driftOnly = phase_err_unwrap * dx;

    % raw shift keeps the base translation plus the wrapped correction
    shift_locked = expected_mod + phase_err_sites;
    shift_locked_unwrap = unwrap_nearest_local(shift_locked, N_wave);
    shift_dist_raw = shift_locked_unwrap * dx;
end

% discard early samples before slope fitting
[dc_raw, se_raw, r2_raw] = fit_slope_local(t_hist, shift_dist_raw);
[dc_do,  se_do,  r2_do ] = fit_slope_local(t_hist, shift_dist_driftOnly);

fprintf('  drift fit: raw shift slope = %.3e, driftOnly slope = %.3e\n', ...
    dc_raw, dc_do);

% collect diagnostics
out = struct();

out.dt          = dt;
out.n_repeats   = n_repeats;
out.M           = M;
out.align_tol   = align_tol;
out.newton_tol  = params.newton_tol;

out.c_ref       = c_ref;
out.tw_res_rel  = tw.tw_res_rel;
out.nsoli_ierr  = tw.nsoli_ierr;
out.nsoli_hist  = tw.nsoli_hist;

out.delta_c_raw       = dc_raw;
out.delta_c_driftOnly = dc_do;
out.delta_c_raw_se    = se_raw;
out.delta_c_raw_r2    = r2_raw;
out.delta_c_se        = se_do;
out.delta_c_r2        = r2_do;

[out.delta_c_first_half_raw, out.delta_c_second_half_raw] = ...
    split_half_fit_local(t_hist, shift_dist_raw);
[out.delta_c_first_half_driftOnly, out.delta_c_second_half_driftOnly] = ...
    split_half_fit_local(t_hist, shift_dist_driftOnly);

[out.corr_shift_Ealign_raw, out.corr_dshift_Ealign_raw] = ...
    corr_metrics_local(t_hist, shift_dist_raw, E_align);
[out.corr_shift_Ealign_driftOnly, out.corr_dshift_Ealign_driftOnly] = ...
    corr_metrics_local(t_hist, shift_dist_driftOnly, E_align);

out.E_align_max    = max(E_align);
out.E_direct_max   = max(E_direct);
out.energy_rel_max = max(Erel_hist);

out.omega_max_eq        = omega_max_eq;
out.dt_omega_eq         = dt * omega_max_eq;
out.dt_omega_eq_over_pi = out.dt_omega_eq / pi;

out.omega_eff0            = omega_eff0;
out.omega_eff_max         = omega_eff_max;
out.dt_omega_eff0         = dt * omega_eff0;
out.dt_omega_eff_max      = dt * omega_eff_max;
out.dt_omega_eff0_over_pi = out.dt_omega_eff0 / pi;
out.dt_omega_eff_max_over_pi = out.dt_omega_eff_max / pi;

out.T_period    = T_period;
out.steps       = steps;
out.store_every = store_every;
out.dt_sample   = dt_sample;

out.t_hist               = t_hist;
out.E_align              = E_align;
out.E_direct             = E_direct;
out.shift_best           = shift_best;
out.shift_dist_raw       = shift_dist_raw;
out.shift_dist_driftOnly = shift_dist_driftOnly;
out.Erel_hist            = Erel_hist;

if MAKE_PLOTS
    figure;
    semilogy(t_hist, E_direct, 'LineWidth', 1.2); hold on;
    semilogy(t_hist, E_align,  'LineWidth', 1.2);
    grid on;
    xlabel('time');
    ylabel('relative error');
    legend({'direct', 'aligned'}, 'Location', 'best');

    figure;
    plot(t_hist, shift_dist_driftOnly, 'LineWidth', 1.2); hold on;
    j0 = max(2, round(0.1*numel(t_hist)));
    p = polyfit(t_hist(j0:end), shift_dist_driftOnly(j0:end), 1);
    plot(t_hist, polyval(p, t_hist), '--', 'LineWidth', 1.2);
    grid on;
    xlabel('time');
    ylabel('shift (distance)');
    legend({'driftOnly', 'linear fit'}, 'Location', 'best');

    figure;
    semilogy(t_hist, Erel_hist, 'LineWidth', 1.2);
    grid on;
    xlabel('time');
    ylabel('relative energy drift');
end

end

function [dc_fit, dc_se, dc_r2] = fit_slope_local(t_hist, y_hist)

dc_fit = NaN;
dc_se  = NaN;
dc_r2  = NaN;

if numel(t_hist) < 20
    return;
end

j0 = max(2, round(0.1 * numel(t_hist)));
tt = t_hist(j0:end);
yy = y_hist(j0:end);

p = polyfit(tt, yy, 1);
dc_fit = p(1);

yhat = polyval(p, tt);
r = yy - yhat;

SSE = sum(r.^2);
SST = sum((yy - mean(yy)).^2);
dc_r2 = 1 - SSE / max(SST, eps);

nfit = numel(tt);
s2 = SSE / max(nfit - 2, 1);
Sxx = sum((tt - mean(tt)).^2);

dc_se = sqrt(s2 / max(Sxx, eps));
end

function [dc1, dc2] = split_half_fit_local(t, y)

dc1 = NaN;
dc2 = NaN;

if numel(t) < 40
    return;
end

j0 = max(2, round(0.1 * numel(t)));
tt = t(j0:end);
yy = y(j0:end);

mid = floor(numel(tt) / 2);
p1 = polyfit(tt(1:mid), yy(1:mid), 1);
p2 = polyfit(tt(mid+1:end), yy(mid+1:end), 1);

dc1 = p1(1);
dc2 = p2(1);
end

function [c1, c2] = corr_metrics_local(t, shift, Ealign)

c1 = NaN;
c2 = NaN;

if numel(t) < 5
    return;
end

C = corrcoef(shift(:), Ealign(:));
if all(isfinite(C(:)))
    c1 = C(1,2);
end

dshift = gradient(shift(:), t(:));
C2 = corrcoef(dshift, Ealign(:));
if all(isfinite(C2(:)))
    c2 = C2(1,2);
end
end

function [xdot, vdot] = rhs_local(x, v, F)

n = length(x);

xdot = zeros(n,1);
vdot = zeros(n,1);

xdot(2:n-1) = v(2:n-1);
vdot(2:n-1) = F(x(3:n)-x(2:n-1)) - F(x(2:n-1)-x(1:n-2));

xdot(1) = v(1);
vdot(1) = F(x(2)-x(1)) - F(x(1)-x(n));

xdot(n) = v(n);
vdot(n) = F(x(1)-x(n)) - F(x(n)-x(n-1));
end

function [best_error, best_delta] = align_error_frac_local(x_final, x_initial, N_wave, q_wave, align_tol)

opts = optimset('TolX', align_tol, 'Display', 'off');

xf = x_final(1:N_wave);
xi = x_initial(1:N_wave);

xf0 = xf - mean(xf);
xi0 = xi - mean(xi);

obj = @(delta) relerr_frac_local(xf0, xi0, delta, q_wave);

half = floor(N_wave / 2);
best_e = inf;
best_int = 0;

for s = -half:half
    e = obj(s);
    if e < best_e
        best_e = e;
        best_int = s;
    end
end

a = max(best_int - 1, -N_wave/2);
b = min(best_int + 1,  N_wave/2);

[best_delta, ~] = fminbnd(obj, a, b, opts);

xs = fracshift_fft_local(xf, best_delta, q_wave);
den = norm(xi);
if den == 0
    den = 1;
end

best_error = norm(xs - xi) / den;
end

function e = relerr_frac_local(xf, xi, delta, q_wave)

xs = fracshift_fft_local(xf, delta, q_wave);
den = norm(xi);
if den == 0
    den = 1;
end

e = norm(xs - xi) / den;
end

function xs = fracshift_fft_local(x, delta, q_wave)

x = x(:);
q_wave = q_wave(:);

X  = fft(x);
Xs = X .* exp(-1i * q_wave * delta);

xs = ifft(Xs, 'symmetric');
xs = xs(:);
end

function s_unwrap = unwrap_nearest_local(s, N_wave)

s = s(:);
s_unwrap = zeros(size(s));

if isempty(s)
    return;
end

s_unwrap(1) = s(1);
for m = 2:numel(s)
    k = round((s_unwrap(m-1) - s(m)) / N_wave);
    s_unwrap(m) = s(m) + k * N_wave;
end
end

function H = energy_FPU_local(x, v, Phi)

dx = x([2:end 1]) - x;
H = 0.5 * sum(v.^2) + sum(Phi(dx));
end

function dp0 = estimate_phip_prime0_local(phip)

h = 1e-6;
dp0 = (phip(h) - phip(-h)) / (2*h);

if ~isfinite(dp0)
    dp0 = 1;
end
end

function y = wrap_to_half_local(x, N)

y = mod(x + N/2, N) - N/2;

tol = 1e-12;
mask = abs(y + N/2) < tol;
y(mask) = +N/2;
end

function kp = phip_prime_vec_local(phip, s)

s = s(:);
h = 1e-6 * max(1, abs(s));

kp = (phip(s + h) - phip(s - h)) ./ (2*h);

bad = ~isfinite(kp);
if any(bad)
    kp(bad) = estimate_phip_prime0_local(phip);
end
end
