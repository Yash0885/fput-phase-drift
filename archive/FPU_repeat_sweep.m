% FPU stability sweep with repeating waves
% doing: continuation -> repeat -> Verlet -> align and compare

close all
clear

% force + potential
phip = @(s) s + s.^3;
Phi  = @(s) 0.5*s.^2 + 0.25*s.^4;

L = 16;
n = [-L:L-1]';
k = pi/L;
N = length(n);
q = k*fftshift([-N/2:N/2-1]');

% run length (periods)
n_periods = 1000;

% how many copies of the wave
n_repeats = 3;

% sh values to test
shift_values = [1e-4 2e-4 5e-4 ...
                1e-3 2e-3 5e-3 ...
                1e-2 1.5e-2 2e-2 2.5e-2 3e-2 4e-2 ...
                5e-2 7.5e-2 1e-1];

% time step
dt = 5e-3;

% storing results
best_error_list = nan(size(shift_values));
best_shift_list = nan(size(shift_values));
c_list          = nan(size(shift_values));
stable_list     = false(size(shift_values));
energy_drift    = nan(size(shift_values));

figure('Position', [100, 100, 1400, 800]);

% reuse continuation guess
uu = 0.25*cos(k*n);
c0 = sqrt(2*(1-cos(k)))/k;
sh_prev = 0;

for idx = 1:length(shift_values)

    sh_target = shift_values(idx);

    % continuing from previous solution
    sh = sh_prev;
    while sh < sh_target
        dsh = min(1e-4 + 0.05*sh, 5e-3);
        sh  = min(sh + dsh, sh_target);

        c = c0 + sh;
        uu = nsoli(uu, @(u)FPU_rhs(u,phip,q,c,k), [1e-12,1e-8], [40,40,0.9,3,20]);
    end
    sh_prev = sh_target;

    c = c0 + sh_target;
    c_list(idx) = c;

    % Initial velocity
    vper = ifft(-c*(1i*q).*fft(uu),'symmetric');

    % repeating the wave
    x = repmat(uu, n_repeats, 1);
    v = repmat(vper, n_repeats, 1);
    x_initial = x;

    % time settings
    T_period = 2*pi/(c*k);
    T = n_periods * T_period;
    steps = round(T/dt);

    % sampling energy about once per period
    store_every = max(1, round(T_period/dt));

    % initial acceleration
    [~, a] = rhs(x, v, phip);

    % initial energy
    E0 = energy_FPU(x, v, Phi);
    Erel_max = 0;

    % Velocity Verlet 
    for step = 1:steps
        x_new = x + v*dt + 0.5*dt^2*a;
        [~, a_new] = rhs(x_new, v, phip);
        v_new = v + 0.5*(a + a_new)*dt;

        x = x_new;
        v = v_new;
        a = a_new;

        % tracking energy drift
        if mod(step, store_every) == 0
            E_now = energy_FPU(x, v, Phi);
            Erel = abs(E_now - E0) / max(1, abs(E0));
            if Erel > Erel_max
                Erel_max = Erel;
            end
        end
    end

    x_final = x;

    % aligning final wave to initial
    N_wave = length(uu);
    [best_error, best_shift, x_aligned] = align_error(x_final, x_initial, N_wave);

    best_error_list(idx) = best_error;
    best_shift_list(idx) = best_shift;
    energy_drift(idx)    = Erel_max;

    % arbitrary stability rule
    is_stable = (best_error < 1e-2) && (Erel_max < 1e-3);
    stable_list(idx) = is_stable;

    fprintf('sh=%.6g | c=%.8f | err=%.3e | shift=%d | energy=%.2e | stable=%d\n', ...
        sh_target, c, best_error, best_shift, Erel_max, is_stable);

    % plotting comparison
    subplot(4,5,idx);
    n_plot = [-(n_repeats*L):(n_repeats*L)-1]';

    plot(n_plot, x_initial, 'b-', 'LineWidth', 2); hold on;
    plot(n_plot, x_aligned, 'r--', 'LineWidth', 1.5);

    if is_stable
        title(sprintf('sh=%.3g stable', sh_target),'Color','g');
    else
        title(sprintf('sh=%.3g unstable', sh_target),'Color','r');
    end

    grid on;
    xlabel('position');
    ylabel('u');

    text(0.5,0.95,sprintf('err=%.2e, shift=%d, E=%.1e',best_error,best_shift,Erel_max),...
        'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',7);

    if idx==1
        legend({'initial','final aligned'});
    end

    drawnow;
end

sgtitle(sprintf('wave stability after %d periods (%d repeats)', n_periods, n_repeats));

% plots
figure;
plot(shift_values, best_error_list,'o-','LineWidth',1.5);
set(gca,'YScale','log'); grid on;
xlabel('sh'); ylabel('aligned error');

figure;
plot(shift_values, energy_drift,'o-','LineWidth',1.5);
set(gca,'YScale','log'); grid on;
xlabel('sh'); ylabel('relative energy drift');

figure;
plot(shift_values, best_shift_list,'o-','LineWidth',1.5);
grid on;
xlabel('sh'); ylabel('best shift');

% functions

function res = FPU_rhs(u, phip, q, c, k)
    ukp1 = ifft(exp(+1i*q).*fft(u),'symmetric');
    ukm1 = ifft(exp(-1i*q).*fft(u),'symmetric');
    upp  = ifft((1i*q).^2.*fft(u),'symmetric');
    res  = -(c)^2*upp + (phip(ukp1-u) - phip(u-ukm1));
end

function [xdot, vdot] = rhs(x, v, F)
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

function [best_error, best_shift, x_best] = align_error(x_final, x_initial, N_wave)
    half = floor(N_wave/2);
    best_error = inf;
    best_shift = 0;
    x_best = x_final;

    den = norm(x_initial);
    if den == 0
        den = 1;
    end
    % Check error of all shifts to find best shift 
    for shift_idx = -half:half
        x_shifted = circshift(x_final, shift_idx);
        err = norm(x_shifted - x_initial) / den;

        if err < best_error
            best_error = err;
            best_shift = shift_idx;
            x_best = x_shifted;
        end
    end
end

function H = energy_FPU(x, v, Phi)
    dx = x([2:end 1]) - x; % x([2:end 1]) means everything is shifted left
    H = 0.5*sum(v.^2) + sum(Phi(dx));
end
