%% Dual Koopman (kernel sections in an RKHS) for ensemble sea-level data
%  Methods compared: DMD, kEDMD, SpecRKHS-Obs (verified dual Koopman eigenpairs)
%
%  Needs in workspace : EnVfull_Xa_36, EnVfull_Ya_36  (n x M*T_snap, Ya(:,k) = state one step after Xa(:,k))
%  Needs on path      : generate_matrices_kernelized, verified_eigenvalues, redblueTecplot, cmap.mat (cmap2)
%  MATLAB             : R2022a+ (clim), R2020a+ (exportgraphics). Use caxis instead of clim on older versions.
%
%  [FIX] tags mark places where the original code was changed.

rng(0)

%% ===================================================================
%  1. DATA AND TRAIN / TEST SPLIT
%  ===================================================================
Xa_all = EnVfull_Xa_36;
Ya_all = EnVfull_Ya_36;
n         = size(Xa_all, 1);            % state dimension (3880)   [FIX] was called p
M         = 10;                         % number of ensemble members
T_snap    = size(Xa_all, 2) / M;        % snapshot pairs per member
grid_size = [97, 40];
%assert(n == prod(grid_size), 'State dimension does not match grid_size.');
nLAND = (1:n).';                        % indices of valid grid points

T_train_split = 29;                     % snapshot pairs used for training (per member)
steps         = T_snap - T_train_split; % forecast horizon                [FIX] T_len removed (duplicate of split)

% [FIX] check that Ya really is Xa shifted by one step inside every member
shift_err = zeros(1, M);
for m = 1:M
    cs = (m-1)*T_snap + 1;
    A_ = Xa_all(:, cs+1 : cs+T_snap-1);
    B_ = Ya_all(:, cs   : cs+T_snap-2);
    shift_err(m) = norm(A_ - B_, 'fro') / norm(B_, 'fro');
end
if max(shift_err) > 1e-8
    warning('Ya is not Xa shifted by one step (max rel. mismatch %.2e). Check the indexing below.', max(shift_err));
end

X_train_list = cell(1, M);  Y_train_list = cell(1, M);
Y_test_list  = cell(1, M);
for m = 1:M
    cs = (m-1)*T_snap + 1;
    X_train_list{m} = Xa_all(:, cs : cs+T_train_split-1);
    Y_train_list{m} = Ya_all(:, cs : cs+T_train_split-1);
    % [FIX] Y_test_list{m}(:,1) is the LAST TRAINING STATE x0 (= Y_train_list{m}(:,end));
    %       columns 2:end are the `steps` future states to be forecast.
    %       The original used Ensembles{m}(:,29), which (if column k = state k) is one step
    %       BEFORE the last training state and puts training data inside the "forecast".
    Y_test_list{m}  = Ya_all(:, cs+T_train_split-1 : cs+T_snap-1);
end
assert(size(Y_test_list{1}, 2) == steps + 1);

X_train = horzcat(X_train_list{:});     % n x N   (N = 29*10 = 290)
Y_train = horzcat(Y_train_list{:});
N       = size(X_train, 2);

%% ===================================================================
%  2. KERNEL MATRICES AND VERIFIED EIGENPAIRS
%  ===================================================================
ker = @kernel_matern;
[G, A, R] = generate_matrices_kernelized(X_train, Y_train, ker);

[Lambda_res, F_res, Lambda, F_all, res, res_verif, idx_verified, W, W_res] = ...
    verified_eigenvalues(G, A, R, N);                       % [FIX] num = N instead of hard-coded 290
fprintf('%d verified eigenvalues (out of %d)\n', length(idx_verified), numel(Lambda));

%% ===================================================================
%  3. PERRON-FROBENIUS MODES
%  ===================================================================
% [FIX] lon/lat: the original used data.y / data.z (from a different dataset) -> undefined here.
lon = 1:grid_size(1);   lat = 1:grid_size(2);     % <-- replace by the real coordinates
                                                   %     (assumes the 97-dimension is longitude; swap if not)
for k_mode = [2 4]                                 % [FIX] loop variable no longer overwrites idx/F/r
    if k_mode > numel(Lambda_res), continue; end
    L_k = Lambda_res(k_mode);
    F_k = F_res(:, k_mode);
    Phi = ((G*F_k) \ X_train.').';                 % n x 1 spatial mode

    v = nan(n, 1);   v(nLAND) = Phi(:);
    v = reshape(v, grid_size);

    figure; ax = axes;                             % [FIX] plain axes instead of nexttile
    imagesc(lon, lat, abs(v).', 'AlphaData', ~isnan(v).');
    colormap(ax, brighten(redblueTecplot(21), -0.55));
    colorbar;
    set(ax, 'Color', [1 1 1]*0.6, 'FontSize', 10);
    axis xy; axis tight;                           % [FIX] removed xlim/ylim [0 2e4] and axis equal (other dataset)
    xlabel('longitude index', 'FontSize', 10); ylabel('latitude index', 'FontSize', 10);
    title(sprintf('$|\\Phi_{%d}|,\\ \\lambda=%.4f%+.4fi$', k_mode, real(L_k), imag(L_k)), ...
          'Interpreter', 'latex', 'FontSize', 12);
end

%% ===================================================================
%  4. SPURIOUS AND VERIFIED EIGENVALUES
%  ===================================================================
figure
scatter(angle(Lambda), log(abs(Lambda)), 200, res, '.', 'LineWidth', 1);
hold on
scatter(angle(Lambda_res), log(abs(Lambda_res)), 500, res_verif, '.', 'LineWidth', 1);
box on
clim([0, 0.1])
load('cmap.mat')
colormap(cmap2); colorbar
xlabel('$\mathrm{arg}(\lambda)$', 'interpreter', 'latex', 'fontsize', 18)
ylabel('$\mathrm{log}(|\lambda|)$', 'interpreter', 'latex', 'fontsize', 18)
title(['Residuals for sea level data', newline], 'interpreter', 'latex', 'fontsize', 18)
ax = gca; ax.FontSize = 18; axis([-pi pi -20e-3 1e-3])
for k = -5:5
    plot(k*pi/6*[1 1], [-20e-3 1e-3], '--', 'Color', 'black')
end
xticks(-pi:pi/6:pi)
ax.TickLabelInterpreter = 'latex';                 % [FIX] set on this axes, not globally via groot
xticklabels({'$-\pi$','','$-2\pi/3$','','$-\pi/3$','','$0$','','$\pi/3$','','$2\pi/3$','','$\pi$'})
xtickangle(30)
exportgraphics(gcf, 'sea_level_evals_angle.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

% (pseudospectrum block from the original is unchanged and left commented out in your file)

%% ===================================================================
%  5. FORECAST MACHINERY (everything that does NOT depend on the ensemble member)
%  ===================================================================
kvec = @(z) ker(X_train, z).';                     % N x 1 vector k(x_i, z); [FIX] replaces the for-loops

g_train = mean(X_train, 1).';                      % observable g(x) = spatial mean, at the N training points

% ---- SpecRKHS-Obs (dual Koopman: expand k_{x0} in the dual eigenfunctions) ----
Lam           = Lambda_res(:);
assert(size(F_res,2) == numel(Lam));
GF            = G * F_res;                         % N x nres
Xmodes_kmd    = F_res.' * X_train.';               % nres x n     [FIX] was F_res.'*X_train (290 vs 3880 -> error)
gmodes_kmd    = F_res.' * g_train;                 % nres x 1     [FIX] observable modes fitted on ALL N training points

% ---- kEDMD ----
use_conj_kedmd = true;                             % convention taken from the original code; see sanity check below
LamK = Lambda(:);   if use_conj_kedmd, LamK = conj(LamK); end
assert(size(W,2) == numel(LamK), 'W and Lambda sizes are inconsistent.');
GW            = G * W;
Xmodes_kedmd  = GW \ X_train.';                    % [FIX] modes fitted on X_train (with exponent 1:steps), not Y_train
gmode_kedmd   = GW \ g_train;                      % [FIX] 290-point fit, not 29/30 points of one member

% ---- standard DMD (computed ONCE, not inside the member loop) ----
[U, S, ~]  = svd(X_train, 'econ');
r_dmd      = rank(S);                              % consider truncating: r_dmd = N makes K nearly exactly-determined
U          = U(:, 1:r_dmd);
PXs        = X_train.' * U;
PYs        = Y_train.' * U;
K_dmd      = PXs \ PYs;
[W_dmd, LAM_dmd] = eig(K_dmd, 'vector');
PXr        = PXs * W_dmd;
Xmodes_dmd = PXr \ X_train.';                      % [FIX] the original 'c' had inconsistent dimensions

% Forecast handles. x0: n x 1 initial state, tt: row vector of lead times, lam: eigenvalues used.
forecast_kmd   = @(x0,tt,lam) real( Xmodes_kmd.'   * ( (GF \ kvec(x0)) .* (lam(:).^tt) ) );
forecast_kedmd = @(x0,tt,lam) real( Xmodes_kedmd.' * ( (W.' * kvec(x0)) .* (lam(:).^tt) ) );
forecast_dmd   = @(x0,tt)     real( Xmodes_dmd.'   * ( (W_dmd.' * (U.' * x0)) .* (LAM_dmd(:).^tt) ) );
obs_kmd        = @(x0,tt,lam) real( ( (GF \ kvec(x0)) .* (lam(:).^tt) ).' * gmodes_kmd ).';
obs_kedmd      = @(x0,tt,lam) real( ( (W.' * kvec(x0)) .* (lam(:).^tt) ).' * gmode_kedmd ).';
kern_kmd       = @(x0,tt,lam) real( GF * ( (GF \ kvec(x0)) .* (lam(:).^tt) ) );   % predicted k_{x_t}(x_i)
rel_err_cols   = @(P,Q) sum(abs(P - Q).^2, 1) ./ sum(abs(Q).^2, 1);

%% ===================================================================
%  6. SANITY CHECKS  (run these before trusting any forecast plot)
%  ===================================================================
nrm = @(a,b) norm(a-b) / norm(b);
xc = X_train(:,1);   yc = Y_train(:,1);
fprintf('\nOne-step check on training pair 1 (relative error, should be << 1):\n');
fprintf('  SpecRKHS-Obs, lambda       : %.3e\n', nrm(forecast_kmd(xc,1,Lam),         yc));
fprintf('  SpecRKHS-Obs, conj(lambda) : %.3e\n', nrm(forecast_kmd(xc,1,conj(Lam)),   yc));
fprintf('  kEDMD,        lambda       : %.3e\n', nrm(forecast_kedmd(xc,1,Lambda(:)),       yc));
fprintf('  kEDMD,        conj(lambda) : %.3e\n', nrm(forecast_kedmd(xc,1,conj(Lambda(:))), yc));
fprintf('  DMD                        : %.3e\n', nrm(forecast_dmd(xc,1),             yc));
fprintf('  Observable == mean of state forecast (linear g): %.3e\n\n', ...
        abs(obs_kmd(xc,1,Lam) - mean(forecast_kmd(xc,1,Lam),1)));
% If the conj(lambda) variant is clearly better for a method, flip the convention for that method.

%% ===================================================================
%  7. FORECASTS AND ERRORS FOR EVERY ENSEMBLE MEMBER
%  ===================================================================
tt = 1:steps;      t0 = 0:steps;

err_dmd_all    = zeros(M, steps);
err_kedmd_all  = zeros(M, steps);
err_kmd_all    = zeros(M, steps);
err_kern_all   = zeros(M, steps+1);
w_exact_all    = zeros(M, steps+1);
w_kmd_all      = zeros(M, steps+1);
w_kedmd_all    = zeros(M, steps+1);

for m = 1:M
    real_window = Y_test_list{m};                  % n x (steps+1): x0, x1, ..., x_steps
    x0          = real_window(:, 1);
    real_data   = real_window(:, 2:end);           % n x steps

    % --- full-state forecasts and relative errors (per lead time) ---
    err_kmd_all(m,:)   = rel_err_cols(forecast_kmd(x0, tt, Lam),   real_data);
    err_kedmd_all(m,:) = rel_err_cols(forecast_kedmd(x0, tt, LamK), real_data);
    err_dmd_all(m,:)   = rel_err_cols(forecast_dmd(x0, tt),        real_data);

    % --- spatial-mean observable ---
    w_exact_all(m,:) = mean(real_window, 1);
    w_kmd_all(m,:)   = obs_kmd(x0, t0, Lam);
    w_kedmd_all(m,:) = obs_kedmd(x0, t0, LamK);

    % --- kernel sections: true k_{x_t}(x_i) vs predicted (K*)^t k_{x0} ---
    K_true = zeros(N, steps+1);
    for t = 1:steps+1
        K_true(:,t) = kvec(real_window(:,t));
    end
    err_kern_all(m,:) = rel_err_cols(kern_kmd(x0, t0, Lam), K_true);
end

%% ===================================================================
%  8. PLOTS (ensemble averaged)
%  ===================================================================
c_dmd = [0 0.4470 0.7410];  c_ked = [0.8500 0.3250 0.0980];  c_kmd = [0.9290 0.6940 0.1250];

% ---- 8a. full-state relative forecast error ----
figure; hold on
plot(tt, mean(err_dmd_all,1),   'linewidth', 2)
plot(tt, mean(err_kedmd_all,1), 'linewidth', 2)
plot(tt, mean(err_kmd_all,1),   'linewidth', 2)
grid on; box on
title('Relative forecast errors comparison (ensemble average)', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Steps)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend({'DMD','kEDMD','SpecRKHS-Obs'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'best')
ax = gca; ax.FontSize = 18;
exportgraphics(gcf, 'sea_level_error_delay.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

% ---- 8b. spatial-mean observable: exact vs kEDMD vs SpecRKHS-Obs ----
weather_exact = mean(w_exact_all, 1);
weather_kmd   = mean(w_kmd_all,   1);
weather_kedmd = mean(w_kedmd_all, 1);

figure; hold on
p1 = plot(t0, weather_kedmd, 'linewidth', 2, 'color', c_ked);
p2 = plot(t0, weather_kmd,   'linewidth', 2, 'color', c_kmd);
p3 = plot(t0, weather_exact, '--', 'linewidth', 2, 'color', c_dmd);
grid on; box on
title('Weather Ensemble Field Forecast Comparison', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Steps)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Spatial Mean Observable Value', 'interpreter', 'latex', 'fontsize', 18)
legend([p3 p1 p2], {'Exact','kEDMD','SpecRKHS-Obs'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'northeast')
ax = gca; ax.FontSize = 18;
exportgraphics(gcf, 'weather_observable_prediction_mean.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

% ---- 8c. same without kEDMD ----
figure; hold on
p2 = plot(t0, weather_kmd,   'linewidth', 2, 'color', c_kmd);
p3 = plot(t0, weather_exact, 'linewidth', 2, 'color', c_dmd);
grid on; box on
title('Weather Ensemble Forecast Comparison', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Steps)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Spatial Mean Weather Value', 'interpreter', 'latex', 'fontsize', 18)
legend([p3 p2], {'Exact','SpecRKHS-Obs'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'southwest')
ax = gca; ax.FontSize = 18;
exportgraphics(gcf, 'weather_prediction_nokedmd.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

% ---- 8d. relative error of the observable forecast ----
% [FIX] The original Figs. 2-4 used `weather_observable_exact`, which was never defined. The observable IS the
%       spatial mean, so "exact observable" and "exact mean weather" are the same curve; the duplicate figures
%       are merged into this one.
% NOTE: if the spatial mean crosses zero, this relative error blows up. In that case divide by max(abs(exact))
%       or by the RMS of `weather_exact` instead.
er_obs = abs(weather_exact - weather_kmd).^2 ./ max(abs(weather_exact).^2, eps);
figure
semilogy(t0, er_obs, 'linewidth', 2, 'color', c_kmd)
grid on
title('Weather Observable Relative Forecast Error', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Steps)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend({'SpecRKHS-Obs'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'northeast')
exportgraphics(gcf, 'weather_error_observable.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

% ---- 8e. relative error of the kernel-section forecast ----
% [FIX] Kx0_future_vals / Kx0_predict were undefined; now K_true / kern_kmd, evaluated at the N training points.
figure
semilogy(t0, mean(err_kern_all, 1), 'linewidth', 2, 'color', c_kmd)
grid on
title('Relative Forecast Error for Kernel Function', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Steps)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend({'SpecRKHS-Obs'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'northeast')
exportgraphics(gcf, 'weather_prediction_kernel_error.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% ===================================================================
%  LOCAL FUNCTIONS (must stay at the end of the script)
%  ===================================================================
% Matern kernel  k(r) = (sigma r)^nu K_nu(sigma r),  nu = 2.   k(0) = 2^(nu-1) Gamma(nu) = 2.
% x: n x m (or n x 1), t: n x 1 (or n x m); returns 1 x m.
% [FIX] besselk(2,.) instead of besselk(-2,.) (identical, K_{-nu} = K_nu, but clearer);
%       the old comment "d=60330" did not match this problem (n = 3880).
function ker = kernel_matern(x, t)
    sigma = 1/10000;                 % length-scale: worth tuning, e.g. sigma = 1/median(pairwise distances)
    r   = vecnorm(x - t);
    ker = zeros(1, size(r, 2));
    pos = r > 0;
    ker(pos)  = (sigma*r(pos)).^2 .* besselk(2, sigma*r(pos));
    ker(~pos) = 2;
end