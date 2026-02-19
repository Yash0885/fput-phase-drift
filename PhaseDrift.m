% Verify phase drift: best shift vs time (sampling once per period)
% The goal is to check if shift grows linearly while shape error stays small

close all
clear

phip = @(s) s + s.^3;

L = 16;
n = [-L:L-1]';
k = pi/L;
N = length(n);
q = k*fftshift([-N/2:N/2-1]');

dt = 5e-3;

n_repeats = 3;

% pick one case to verify
sh_target = 2e-2;       % try: 1e-4, 1e-2, 2e-2, 3e-2, 5e-2
n_periods = 1000;

% traveling wave by continuation 
uu = 0.25*cos(k*n);

c0 = sqrt(2*(1-cos(k)))/k;
sh = 0;

while sh < sh_target
    dsh = min(1e-4 + 0.05*sh, 5e-3);
    sh  = min(sh + dsh, sh_target);

    c = c0 + sh;
    uu = nsoli(uu, @(u)FPU_rhs(u,phip,q,c,k), [1e-12,1e-8], [40, 40, 0.9, 3, 20]);
end

c = c0 + sh_target;

vper = ifft(-c*(1i*q).*fft(uu), 'symmetric');

if max(abs(uu)) < 1e-6
    error('Trivial uu (near zero).')
end

% IC
x = repmat(uu, n_repeats, 1);
v = repmat(vper, n_repeats, 1);
x0 = x;

% integration setup
T_period = 2*pi/(c*k);
T = n_periods * T_period;
steps = round(T/dt);

store_every = max(1, round(T_period/dt));   % once per period
nstore = floor(steps/store_every) + 1;

t_hist   = zeros(nstore,1);
err_hist = zeros(nstore,1);
shf_hist = zeros(nstore,1);

store_i = 1;
t_hist(store_i) = 0;
[err_hist(store_i), shf_hist(store_i)] = align_error(x0, x0, length(uu));

[~, a] = rhs(x, v, phip);

% Verlet + sampling
for step = 1:steps
    x_new = x + v*dt + 0.5*dt^2*a;
    [~, a_new] = rhs(x_new, v, phip);
    v_new = v + 0.5*(a + a_new)*dt;

    x = x_new;
    v = v_new;
    a = a_new;

    if mod(step, store_every) == 0
        store_i = store_i + 1;
        t_hist(store_i) = step*dt;
        [err_hist(store_i), shf_hist(store_i)] = align_error(x, x0, length(uu));
    end
end

t_hist = t_hist(1:store_i);
err_hist = err_hist(1:store_i);
shf_hist = shf_hist(1:store_i);

% unwrap shifts so drift can exceed +-N_wave/2
N_wave = length(uu);
shf_unwrap = shf_hist;

for m = 2:numel(shf_unwrap)
    ds = shf_unwrap(m) - shf_unwrap(m-1);
    if ds >  N_wave/2
        shf_unwrap(m:end) = shf_unwrap(m:end) - N_wave;
    elseif ds < -N_wave/2
        shf_unwrap(m:end) = shf_unwrap(m:end) + N_wave;
    end
end

% linear fit for drift rate (ignore early transient)
dc_fit = NaN;
if numel(t_hist) >= 20
    j0 = max(2, round(0.1*numel(t_hist)));
    p = polyfit(t_hist(j0:end), shf_unwrap(j0:end), 1);
    dc_fit = p(1);  % sites / time
end

fprintf('sh=%.4g | c=%.8f | drift rate ~ %.3e (sites/time)\n', sh_target, c, dc_fit);

% plots
figure;
plot(t_hist, shf_unwrap, 'LineWidth', 1.5);
grid on;
xlabel('time t');
ylabel('unwrapped best shift (sites)');
title(sprintf('Phase drift: sh=%.3g, c=%.5f, slope~%.2e', sh_target, c, dc_fit));

figure;
plot(t_hist, err_hist, 'LineWidth', 1.5);
set(gca,'YScale','log');
grid on;
xlabel('time t');
ylabel('best relative error (aligned)');
title('Shape error vs time');

% functions
function res = FPU_rhs(u, phip, q, c, k)
    ukp1 = ifft(exp(+1i*q).*fft(u), 'symmetric');
    ukm1 = ifft(exp(-1i*q).*fft(u), 'symmetric');
    upp  = ifft((1i*q).^2.*fft(u), 'symmetric');
    res  = -(c)^2*upp + (phip(ukp1-u) - phip(u-ukm1));
end

function [xdot, vdot] = rhs(x, v, F)
    n = length(x);
    xdot = zeros(n, 1);
    vdot = zeros(n, 1);

    xdot(2:n-1) = v(2:n-1);
    vdot(2:n-1) = F(x(3:n) - x(2:n-1)) - F(x(2:n-1) - x(1:n-2));

    xdot(1) = v(1);
    vdot(1) = F(x(2) - x(1)) - F(x(1) - x(n));

    xdot(n) = v(n);
    vdot(n) = F(x(1) - x(n)) - F(x(n) - x(n-1));
end

function [best_error, best_shift, x_best] = align_error(x_final, x_initial, N_wave)
    half = floor(N_wave/2);
    best_error = inf;
    best_shift = 0;
    x_best = x_final;

    for shift_idx = -half:half
        x_shifted = circshift(x_final, shift_idx);
        err = norm(x_shifted - x_initial) / norm(x_initial);

        if err < best_error
            best_error = err;
            best_shift = shift_idx;
            x_best = x_shifted;
        end
    end
end
