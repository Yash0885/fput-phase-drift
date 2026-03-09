function tw = build_traveling_wave_frac(params)
% build_traveling_wave_frac
%
% Build one traveling wave at the target speed shift.
% Return the profile, matching velocity, and a residual check.

sh_target    = params.sh_target;
L            = params.L;
k            = params.k;
phip         = params.phip;
newton_tol   = params.newton_tol;
newton_parms = params.newton_parms;

% lattice grid
n = (-L:L-1)';
N = length(n);

% Fourier wave numbers for spectral differentiation
q = k * fftshift((-N/2):(N/2-1))';
q = q(:);

% initial guess for continuation
uu = 0.25*cos(k*n);
uu = uu(:);

% linear reference speed
c0 = sqrt(2*(1 - cos(k))) / k;
sh = 0;

it_hist = [];
ierr = NaN;

% continue in speed shift up to sh_target
while sh < sh_target
    dsh = min(1e-4 + 0.05*sh, 5e-3);
    sh  = min(sh + dsh, sh_target);

    c = c0 + sh;

    % reuse the previous profile as the next Newton guess
    try
        [uu, it_hist, ierr] = nsoli(uu, @(u) FPU_rhs_local(u, phip, q, c), ...
            newton_tol, newton_parms);
    catch
        % some nsoli variants return only the solution
        uu = nsoli(uu, @(u) FPU_rhs_local(u, phip, q, c), ...
            newton_tol, newton_parms);
        it_hist = [];
        ierr = NaN;
    end

    uu = uu(:);

    if ~isnan(ierr) && ierr ~= 0
        warning('nsoli returned ierr=%d at sh=%.6g', ierr, sh);
    end
end

% target speed reused in the time-stepping runs
c_ref = c0 + sh_target;

if max(abs(uu)) < 1e-6
    error('Trivial uu (near zero).');
end

% residual at the target speed
tw_res = FPU_rhs_local(uu, phip, q, c_ref);
tw_res_rel = norm(tw_res) / max(norm(uu), 1);

% consistent initial velocity from u_t ~ -c u_x
vper = ifft(-c_ref * (1i*q) .* fft(uu), 'symmetric');
vper = vper(:);

fprintf('  traveling wave built: max|u|=%.3e, c_ref=%.12g, relres=%.3e\n', ...
    max(abs(uu)), c_ref, tw_res_rel);

% keep traveling-wave data together for the validation runs
tw = struct();
tw.uu         = uu;
tw.vper       = vper;
tw.c_ref      = c_ref;
tw.c0         = c0;
tw.q          = q;
tw.tw_res_rel = tw_res_rel;
tw.nsoli_ierr = ierr;
tw.nsoli_hist = it_hist;

end

function res = FPU_rhs_local(u, phip, q, c)

u = u(:);
q = q(:);

ukp1 = ifft(exp(+1i*q) .* fft(u), 'symmetric');
ukm1 = ifft(exp(-1i*q) .* fft(u), 'symmetric');
upp  = ifft((1i*q).^2 .* fft(u), 'symmetric');

res = -(c)^2 * upp + (phip(ukp1 - u) - phip(u - ukm1));
res = res(:);

end
