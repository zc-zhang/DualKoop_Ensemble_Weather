%plot_sp_dual_kernel_weather.m

% Plot_kernel_ResDMD_weather.m
% Dual Koopman Mode Decomposition for Weather Data
% With sparsity-promoting amplitude selection (NaN-safe)
% Date: 2026/02/05    By Z. Zhang

%clear; close all;

%% ============================================================
% 1. Load and prepare data
%% ============================================================
data_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Ensemble_MSLP';
M = 100;
Xa_ens = cell(M, 1);
Ya_ens = cell(M, 1);

fprintf('Loading %d ensemble members (enMSLP000-099)...\n', M);
for m = 0:M-1
    fname  = sprintf('enMSLP%03d.mat', m);
    S_load = load(fullfile(data_dir, fname));
    Xraw   = S_load.vectorized_MSLP;
    Xa_ens{m+1} = Xraw(:, 2:end-1);  % 20860 × 71
    Ya_ens{m+1} = Xraw(:, 3:end);    % 20860 × 71
end

fprintf('Loading mean trajectory (enMSLP100)...\n');
S_mean    = load(fullfile(data_dir, 'enMSLP100.mat'));
Xraw_mean = S_mean.vectorized_MSLP;
Xa_mean   = Xraw_mean(:, 2:end-1);  % 20860 × 71
Ya_mean   = Xraw_mean(:, 3:end);    % 20860 × 71

%% --- Deviation from ensemble mean ---
Xa = Xa_ens{1,1} - Xa_mean;   % 20860 × 71
Ya = Ya_ens{1,1} - Ya_mean;   % 20860 × 71

%% --- Grid and coastline ---
lat_sub = weatherDat2021AUG_ensemble0.dat_lat;  % 140 × 149
lon_sub = weatherDat2021AUG_ensemble0.dat_lon;  % 140 × 149
latN = 140;
lonN = 149;
shapefile_path = 'D:\Susuki Lab\Testing_Code\data-weather\Data_250401\ne_10m_coastline\ne_10m_coastline.shp';
S = shaperead(shapefile_path);

%% ============================================================
% 2. Kernel EDMD — Dual Koopman
%% ============================================================
%N_modes = 70;

[G, K_star, L, PX, PY, PSI_x] = kernel_ResDMD(...
    Xa, Ya,           ...
    'type', 'Gaussian', ...
    'N',    300,  ...
    'Xb',   Xa        ...
);

[V, D]  = eig(K_star);
Lambda  = diag(D);
nModes  = numel(Lambda);
fprintf('Computed %d dual Koopman modes\n', nModes);

%% ============================================================
% 3. Pre-compute Physical Space Modes (once — outside all loops)
%% ============================================================
fprintf('Computing physical space modes...\n');
Modes_space = zeros(size(Xa,1), nModes);   % 20860 × nModes
for k = 1:nModes
    Modes_space(:,k) = Xa * (PX * V(:,k));
end
fprintf('Modes_space: [%d x %d]\n', size(Modes_space,1), size(Modes_space,2));

%% ============================================================
% 4. Compute Amplitudes
%% ============================================================
fprintf('\n--- Computing dual Koopman amplitudes ---\n');
nTime = size(Xa, 2);   % 71

%% --- Vandermonde matrix (nModes × nTime) ---
Vand = zeros(nModes, nTime);
for k = 1:nModes
    Vand(k,:) = Lambda(k).^(0:nTime-1);
end
fprintf('Vandermonde: [%d x %d]\n', size(Vand,1), size(Vand,2));

%% --- Gram matrix P and projection vector q ---
P_raw = (Modes_space' * Modes_space) .* conj(Vand * Vand');
q_raw = conj(diag(Vand * Xa' * Modes_space));

%% --- Diagnose P ---
fprintf('\n--- P matrix diagnostics ---\n');
fprintf('P is real?        %d\n',   isreal(P_raw));
fprintf('P symmetry error: %.2e\n', norm(P_raw - P_raw', 'fro'));
fprintf('min eigenvalue:   %.4e\n', min(real(eig(P_raw))));
fprintf('condition number: %.4e\n', cond(P_raw));

%% --- Force P to be real symmetric positive definite ---
P_fix = real((P_raw + P_raw') / 2);

% Adaptive regularization scaled to diagonal of P
reg   = max(1e-6 * max(abs(diag(P_fix))), 1e-10);
P_fix = P_fix + reg * eye(size(P_fix));
fprintf('Regularization applied: %.4e\n', reg);
fprintf('min eigenvalue after fix: %.4e\n', min(eig(P_fix)));

% Force q to be real
q_fix = real(q_raw);

%% --- Safe Cholesky ---
try
    Pl = chol(P_fix);
    fprintf('Cholesky succeeded\n');
catch
    fprintf('Cholesky failed — using eig-based fix\n');
    [Qe, De] = eig(P_fix);
    De_vec   = max(real(diag(De)), reg);
    P_fix    = real(Qe * diag(De_vec) * Qe');
    P_fix    = (P_fix + P_fix') / 2;
    Pl       = chol(P_fix);
    fprintf('Cholesky succeeded after eig fix\n');
end

%% --- Least squares amplitudes (dense baseline) ---
b_ls = P_fix \ q_fix;
fprintf('\nLS amplitudes — range: [%.4e, %.4e]\n', ...
    min(abs(b_ls)), max(abs(b_ls)));

%% ============================================================
% 5. Sparsity-Promoting Selection (ADMM)
%% ============================================================
fprintf('\n--- Sparsity-promoting sweep ---\n');

gamma_vals = logspace(-2, 2, 50);
n_active   = zeros(size(gamma_vals));
residuals  = zeros(size(gamma_vals));

for gi = 1:numel(gamma_vals)
    b_sp = sparsifyDMD(P_fix, q_fix, Pl, gamma_vals(gi), nModes);

    if any(isnan(b_sp))
        fprintf('  Warning: NaN at gamma=%.2e — skipping\n', gamma_vals(gi));
        n_active(gi)  = NaN;
        residuals(gi) = NaN;
        continue;
    end

    thr           = 1e-6 * max(abs(b_sp));
    n_active(gi)  = sum(abs(b_sp) > thr);
    residuals(gi) = norm(Xa - Modes_space * diag(b_sp) * Vand, 'fro')^2;

    if mod(gi, 10) == 0
        fprintf('  gamma=%.2e -> %d active modes, residual=%.4e\n', ...
            gamma_vals(gi), n_active(gi), residuals(gi));
    end
end

%% --- Remove NaN entries ---
valid          = isfinite(n_active) & isfinite(residuals);
gamma_valid    = gamma_vals(valid);
n_active_valid = n_active(valid);
resid_valid    = residuals(valid);
fprintf('Valid gamma points: %d / %d\n', sum(valid), numel(gamma_vals));

%% --- Elbow detection ---
r_norm = resid_valid  / max(resid_valid);
n_norm = n_active_valid / max(n_active_valid);
[~, elbow_idx] = min(r_norm + n_norm);
gamma_opt      = gamma_valid(elbow_idx);
fprintf('Optimal gamma* = %.4e -> %d active modes\n', ...
    gamma_opt, n_active_valid(elbow_idx));

%% --- Final sparse amplitudes ---
b_sparse = sparsifyDMD(P_fix, q_fix, Pl, gamma_opt, nModes);

if any(isnan(b_sparse))
    warning('NaN in final b_sparse — falling back to LS amplitudes');
    b_sparse = b_ls;
end

fprintf('Sparse amplitudes — range: [%.4e, %.4e]\n', ...
    min(abs(b_sparse)), max(abs(b_sparse)));

%% --- Sort active modes by sparse amplitude (largest first) ---
active_mask  = abs(b_sparse) > 1e-6 * max(abs(b_sparse));
active_modes = find(active_mask);
[~, amp_order] = sort(abs(b_sparse(active_modes)), 'descend');
eig_list       = active_modes(amp_order);
fprintf('Total active modes: %d\n', numel(eig_list));

%% --- Print top modes summary ---
fprintf('\n%-6s %-10s %-14s %-12s %-12s\n', ...
    'Rank','ModeIdx','|b_sparse|','|Lambda|','Period(dt)');
for r = 1:min(15, numel(eig_list))
    k   = eig_list(r);
    ang = abs(angle(Lambda(k)));
    T   = 2*pi / max(ang, 1e-10);
    fprintf('%-6d %-10d %-14.4e %-12.4f %-12.4f\n', ...
        r, k, abs(b_sparse(k)), abs(Lambda(k)), T);
end

%% ============================================================
% 6. Plot: Sparsity Sweep Elbow Curve
%% ============================================================
figure('Color','w','Position',[100 100 900 400]);
subplot(1,2,1);
semilogx(gamma_valid, n_active_valid, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5);
hold on;
xline(gamma_opt, 'r--', 'LineWidth', 1.5, ...
    'Label', sprintf('\\gamma^*=%.2e', gamma_opt), ...
    'LabelVerticalAlignment','bottom');
hold off;
xlabel('$\gamma$',            'Interpreter','latex', 'FontSize', 12);
ylabel('Active modes',         'FontSize', 12);
title('Sparsity vs $\gamma$', 'Interpreter','latex', 'FontSize', 12);
grid on;

subplot(1,2,2);
semilogx(gamma_valid, resid_valid, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 5);
hold on;
xline(gamma_opt, 'r--', 'LineWidth', 1.5, ...
    'Label', sprintf('\\gamma^*=%.2e', gamma_opt), ...
    'LabelVerticalAlignment','bottom');
hold off;
xlabel('$\gamma$',                              'Interpreter','latex', 'FontSize', 12);
ylabel('$\|X - \hat{X}\|_F^2$',                'Interpreter','latex', 'FontSize', 12);
title('Reconstruction Error vs $\gamma$',       'Interpreter','latex', 'FontSize', 12);
grid on;
sgtitle('Sparsity-Promoting Mode Selection', 'FontSize', 13);

%% ============================================================
% 7. Eigenvalue Spectrum (color = sparse amplitude)
%% ============================================================
figure('Color','w','Position',[100 100 800 600]);
scatter(real(Lambda), imag(Lambda), 60, abs(b_sparse), 'filled', ...
    'MarkerEdgeColor','none');
hold on;
scatter(real(Lambda(eig_list)), imag(Lambda(eig_list)), 120, ...
    abs(b_sparse(eig_list)), 'filled', 'MarkerEdgeColor','r', 'LineWidth', 2);
theta = linspace(0, 2*pi, 200);
plot(cos(theta), sin(theta), 'r--', 'LineWidth', 1.5);
hold off;
colormap(brighten(redblueTecplot(21), -0.55));
cb = colorbar;
ylabel(cb, '$|b_{sparse}|$', 'Interpreter','latex', 'FontSize', 11);
xlabel('Re($\lambda$)', 'Interpreter','latex', 'FontSize', 12);
ylabel('Im($\lambda$)', 'Interpreter','latex', 'FontSize', 12);
title('Dual Koopman Eigenvalues — color = $|b_{sparse}|$', ...
    'Interpreter','latex', 'FontSize', 13);
legend('All modes','Selected modes','Unit circle','Location','best');
grid on; axis equal;

%% ============================================================
% 8. Top 9 Modes — Magnitude |φ|
%% ============================================================
%nPlot    = min(9, numel(eig_list));
nPlot    = numel(eig_list);
%nPlot    = ;
eig_plot = eig_list(1:2:nPlot);

figure('Color','w','Position',[100 100 1200 900]);
tl = tiledlayout(3, 3, 'Padding','compact', 'TileSpacing','compact');
title(tl, 'Top Dual Koopman Modes $|\varphi|$ — Sorted by Sparse Amplitude', ...
    'Interpreter','latex', 'FontSize', 13);

for k = 1:nPlot
    eig_id_i  = eig_plot(k);
    phi_field = abs(reshape(Modes_space(:, eig_id_i), [latN, lonN]));
    ang       = abs(angle(Lambda(eig_id_i)));
    period    = 2*pi / max(ang, 1e-10);

    nexttile;
    p = pcolor(lon_sub, lat_sub, phi_field);
    p.EdgeColor = 'none';
    shading interp;
  colormap(brighten(redblueTecplot(21), -0.55));
cb = colorbar;
    caxis([0, max(phi_field(:))]);
    colorbar;
    set(gca,'YDir','normal');
    axis equal tight;
    hold on;
    for kk = 1:length(S)
        plot(S(kk).X, S(kk).Y, 'k-', 'LineWidth', 0.8);
    end
    hold off;
    xlim([min(lon_sub(:)), max(lon_sub(:))]);
    ylim([min(lat_sub(:)), max(lat_sub(:))]);
    xlabel('Lon (°E)', 'Interpreter','tex', 'FontSize', 9);
    ylabel('Lat (°N)', 'Interpreter','tex', 'FontSize', 9);
    title(sprintf('$V_{%d},\\;|b|=%.2e,\\;|\\lambda|=%.3f,\\;T=%.2f$', ...
        eig_id_i, abs(b_sparse(eig_id_i)), abs(Lambda(eig_id_i)), period), ...
        'Interpreter','latex', 'FontSize', 9);
end

%% ============================================================
% 9. Top 9 Modes — Real Part Re(φ)
%% ============================================================
figure('Color','w','Position',[100 100 1200 900]);
tl2 = tiledlayout(3, 3, 'Padding','compact', 'TileSpacing','compact');
title(tl2, 'Top Dual Koopman Modes Re($\varphi$) — Sorted by Sparse Amplitude', ...
    'Interpreter','latex', 'FontSize', 13);

for k = 1:nPlot
    eig_id_i       = eig_plot(k);
    phi_field_real = real(reshape(Modes_space(:, eig_id_i), [latN, lonN]));
    absMax         = max(abs(phi_field_real(:)));
    ang            = abs(angle(Lambda(eig_id_i)));
    period         = 2*pi / max(ang, 1e-10);

    nexttile;
    p = pcolor(lon_sub, lat_sub, phi_field_real);
    p.EdgeColor = 'none';
    shading interp;
    colormap(gca, brighten(redblueTecplot(21), -0.55));
    caxis([-absMax, absMax]);
    colorbar;
    set(gca,'YDir','normal');
    axis equal tight;
    hold on;
    for kk = 1:length(S)
        plot(S(kk).X, S(kk).Y, 'k-', 'LineWidth', 0.8);
    end
    hold off;
    xlim([min(lon_sub(:)), max(lon_sub(:))]);
    ylim([min(lat_sub(:)), max(lat_sub(:))]);
    xlabel('Lon (°E)', 'Interpreter','tex', 'FontSize', 9);
    ylabel('Lat (°N)', 'Interpreter','tex', 'FontSize', 9);
    title(sprintf('$V_{%d},\\;|b|=%.2e,\\;|\\lambda|=%.3f,\\;T=%.2f$', ...
        eig_id_i, abs(b_sparse(eig_id_i)), abs(Lambda(eig_id_i)), period), ...
        'Interpreter','latex', 'FontSize', 9);
end

%% ============================================================
% 10. Individual 4-panel detail for each top mode
%% ============================================================
for k = 1:nPlot
    eig_id_i = eig_plot(k);
    ang      = abs(angle(Lambda(eig_id_i)));
    period   = 2*pi / max(ang, 1e-10);

    figure('Color','w','Position',[100 100 1400 900]);

    %% (a) Eigenvalue spectrum
    subplot(2,2,1);
    scatter(real(Lambda), imag(Lambda), 60, abs(b_sparse), 'filled', ...
        'MarkerEdgeColor','none');
    hold on;
    scatter(real(Lambda(eig_id_i)), imag(Lambda(eig_id_i)), 150, ...
        'r', 'filled', 'MarkerEdgeColor','k', 'LineWidth', 2);
    theta = linspace(0, 2*pi, 200);
    plot(cos(theta), sin(theta), 'b--', 'LineWidth', 1.5);
    hold off;
 colormap(brighten(redblueTecplot(21), -0.55));
cb = colorbar;
    xlabel('Re($\lambda$)', 'Interpreter','latex', 'FontSize', 11);
    ylabel('Im($\lambda$)', 'Interpreter','latex', 'FontSize', 11);
    title(sprintf('Spectrum — Mode %d  |b|=%.2e', ...
        eig_id_i, abs(b_sparse(eig_id_i))), 'FontSize', 11);
    legend('All (color=|b|)', sprintf('Mode %d', eig_id_i), ...
        'Unit circle', 'Location','best', 'FontSize', 8);
    grid on; axis equal;

    %% (b) RKHS coefficients
    subplot(2,2,2);
    coeffs = abs(V(:, eig_id_i));
    stem(1:length(coeffs), coeffs, 'filled', 'MarkerSize', 4, ...
        'Color',[0.3 0.3 0.8], 'MarkerFaceColor',[0.3 0.3 0.8]);
    hold on;
    stem(eig_id_i, coeffs(eig_id_i), 'filled', 'MarkerSize', 10, ...
        'Color','r', 'MarkerFaceColor','r', 'LineWidth', 2.5);
    xline(eig_id_i, 'r--', 'LineWidth', 1.5, ...
        'Label', sprintf('Mode %d', eig_id_i));
    hold off;
    xlabel('Snapshot Index', 'FontSize', 11);
    ylabel('|Coefficient|',  'FontSize', 11);
    title(sprintf('RKHS Coefficients — Mode %d  |b|=%.2e', ...
        eig_id_i, abs(b_sparse(eig_id_i))), 'FontSize', 11);
    grid on;

    %% (c) Temporal evolution
    subplot(2,2,3);
    t_vec    = (0:nTime-1)';
    phi_j_x0 = PSI_x(1,:) * V(:, eig_id_i);
    temporal  = phi_j_x0 * Lambda(eig_id_i).^t_vec;
    plot(t_vec, real(temporal), 'b-.o', 'LineWidth', 2, 'MarkerSize', 4);
    xlabel('Time step $t$', 'Interpreter','latex', 'FontSize', 11);
    ylabel('Re$(\lambda_j^t \phi_j(x_0))$', 'Interpreter','latex', 'FontSize', 11);
    title(sprintf('$|\\lambda|=%.3f,\\;T=%.2f\\,\\mathrm{dt},\\;|b|=%.2e$', ...
        abs(Lambda(eig_id_i)), period, abs(b_sparse(eig_id_i))), ...
        'Interpreter','latex', 'FontSize', 11);
    grid on;

    %% (d) Spatial mode magnitude
    subplot(2,2,4);
    phi_field = abs(reshape(Modes_space(:, eig_id_i), [latN, lonN]));
    p = pcolor(lon_sub, lat_sub, phi_field);
    p.EdgeColor = 'none';
    shading interp;
colormap(brighten(redblueTecplot(21), -0.55));
cb = colorbar;
    caxis([0, max(phi_field(:))]);
    cb = colorbar;
    ylabel(cb, '$|\varphi|$', 'Interpreter','latex', 'FontSize', 10);
    set(gca,'YDir','normal');
    axis equal tight;
    hold on;
    for kk = 1:length(S)
        plot(S(kk).X, S(kk).Y, 'k-', 'LineWidth', 1.2);
    end
    hold off;
    xlabel('Longitude (°E)', 'Interpreter','tex', 'FontSize', 11);
    ylabel('Latitude (°N)',  'Interpreter','tex', 'FontSize', 11);
    title(sprintf('$V_{%d},\\;|\\lambda|=%.3f,\\;T=%.2f,\\;|b|=%.2e$', ...
        eig_id_i, abs(Lambda(eig_id_i)), period, abs(b_sparse(eig_id_i))), ...
        'Interpreter','latex', 'FontSize', 11);
    xlim([min(lon_sub(:)), max(lon_sub(:))]);
    ylim([min(lat_sub(:)), max(lat_sub(:))]);

    sgtitle(sprintf('Dual Koopman Mode %d — Rank %d by Sparse Amplitude', ...
        eig_id_i, k), 'FontSize', 13, 'FontWeight','bold');
end

%% ============================================================
% 11. 3D Ensemble Trajectories
%% ============================================================
figure('Position',[100 100 1200 800]);
hold on;
n_snapshots      = size(Xa_ens{1}, 2);
snapshots        = 0:n_snapshots-1;
spatial_mean_100 = mean(Xa_mean, 1);
plot3(ones(size(snapshots))*100, snapshots, spatial_mean_100, ...
    'r-', 'LineWidth', 4, 'DisplayName', 'enMSLP100 (Ensemble Mean)');

all_means = zeros(M, n_snapshots);
for m = 1:M
    all_means(m,:) = mean(Xa_ens{m}, 1);
end
cmin = min(all_means(:));
cmax = max(all_means(:));

for m = 1:M
    ensemble_idx = ones(size(snapshots)) * (m-1);
    patch([ensemble_idx, nan], [snapshots, nan], [all_means(m,:), nan], ...
          [all_means(m,:), nan], ...
          'EdgeColor','interp', 'FaceColor','none', 'LineWidth', 1.5);
end
hold off;
xlabel('Ensemble Member Index (0-100)');
ylabel('Snapshot / Time Step (0-70)');
zlabel('Spatial Mean of MSLP');
title('3D Ensemble Trajectories');
colormap('jet'); colorbar;
caxis([cmin cmax]);
legend('Location','best');
grid on;
view(45, 30);
rotate3d on;

%% ============================================================
% 12. Summary
%% ============================================================
fprintf('\n========== Dual Koopman Analysis Summary ==========\n');
fprintf('Data:           %d spatial pts x %d time steps\n', size(Xa,1), size(Xa,2));
fprintf('Modes computed: %d\n', nModes);
fprintf('Reg applied:    %.4e\n', reg);
fprintf('Sparse modes:   %d  (gamma* = %.2e)\n', numel(eig_list), gamma_opt);
fprintf('\nTop %d modes by sparse amplitude:\n', nPlot);
fprintf('%-6s %-10s %-14s %-12s %-12s\n', ...
    'Rank','ModeIdx','|b_sparse|','|Lambda|','Period(dt)');
for r = 1:nPlot
    k   = eig_plot(r);
    ang = abs(angle(Lambda(k)));
    T   = 2*pi / max(ang, 1e-10);
    fprintf('%-6d %-10d %-14.4e %-12.4f %-12.4f\n', ...
        r, k, abs(b_sparse(k)), abs(Lambda(k)), T);
end
fprintf('====================================================\n');




%%%%%%%%%%%
%% === COMPUTE NORMALIZED PERFORMANCE LOSS ===

% % Baseline 1: J_0 = ||X||_F^2  (trivial predictor)
% J_0_data = norm(Xa, 'fro')^2;
% 
% % Baseline 2: J_0 = reconstruction with ALL modes (LS)
% Xa_recon_ls = Modes_space * diag(b_ls) * Vand;
% J_0_ls      = norm(Xa - Xa_recon_ls, 'fro')^2;
% 
% fprintf('||X||_F^2          = %.4e\n', J_0_data);
% fprintf('J(b_ls) full recon = %.4e\n', J_0_ls);
% fprintf('Relative LS error  = %.4f%%\n', sqrt(J_0_ls/J_0_data)*100);
% 
% %% === SWEEP GAMMA with normalized loss ===
% %gamma_vals    = logspace(-5, 5, 50);
% n_active      = zeros(size(gamma_vals));
% residuals_raw = zeros(size(gamma_vals));   % ||X - Xhat||_F^2
% J_gamma       = zeros(size(gamma_vals));   % full loss including L1
% Pi_data       = zeros(size(gamma_vals));   % normalized by ||X||_F^2
% Pi_ls         = zeros(size(gamma_vals));   % normalized by J(b_ls)
% 
% for gi = 1:numel(gamma_vals)
%     gamma = gamma_vals(gi);
%     b_sp  = sparsifyDMD(P_fix, q_fix, Pl, gamma, nModes);
% 
%     if any(isnan(b_sp))
%         residuals_raw(gi) = NaN;
%         J_gamma(gi)       = NaN;
%         Pi_data(gi)       = NaN;
%         Pi_ls(gi)         = NaN;
%         n_active(gi)      = NaN;
%         continue;
%     end
% 
%     % Reconstruction error
%     recon_err         = norm(Xa - Modes_space * diag(b_sp) * Vand, 'fro')^2;
%     residuals_raw(gi) = recon_err;
% 
%     % Full performance loss J_gamma(b) = recon_err + gamma*||b||_1
%     J_gamma(gi)  = recon_err + gamma * sum(abs(b_sp));
% 
%     % Normalized loss — Jovanovic style
%     Pi_data(gi)  = recon_err / J_0_data;    % vs data norm
%     Pi_ls(gi)    = recon_err / max(J_0_ls, 1e-10);  % vs LS solution
% 
%     thr          = 1e-6 * max(abs(b_sp));
%     n_active(gi) = sum(abs(b_sp) > thr);
% end
% 
% %% === FILTER VALID ===
% valid          = isfinite(Pi_data) & isfinite(n_active);
% gamma_valid    = gamma_vals(valid);
% Pi_data_valid  = Pi_data(valid);
% Pi_ls_valid    = Pi_ls(valid);
% n_active_valid = n_active(valid);
% resid_valid    = residuals_raw(valid);
% 
% %% === ELBOW ON NORMALIZED LOSS ===
% % Use Pi_ls — normalized by LS baseline (Jovanovic style)
% r_norm = Pi_ls_valid / max(Pi_ls_valid);
% n_norm = n_active_valid / max(n_active_valid);
% [~, elbow_idx] = min(r_norm + n_norm);
% gamma_opt      = gamma_valid(elbow_idx);
% 
% fprintf('\nOptimal gamma* = %.4e\n', gamma_opt);
% fprintf('Active modes   = %d\n',    n_active_valid(elbow_idx));
% fprintf('Pi(gamma*)     = %.4f%%\n', Pi_ls_valid(elbow_idx)*100);
% 
% %% === PLOT NORMALIZED PERFORMANCE LOSS ===
% figure('Color','w','Position',[100 100 1200 450]);
% 
% subplot(1,3,1);
% semilogx(gamma_valid, n_active_valid, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5);
% hold on;
% xline(gamma_opt, 'r--', 'LineWidth', 1.5, ...
%     'Label', sprintf('\\gamma^*=%.2e', gamma_opt), ...
%     'LabelVerticalAlignment','bottom');
% hold off;
% xlabel('$\gamma$',       'Interpreter','latex', 'FontSize', 12);
% ylabel('Active modes',    'FontSize', 12);
% title('$N_{active}$ vs $\gamma$', 'Interpreter','latex', 'FontSize', 12);
% grid on;
% 
% subplot(1,3,2);
% semilogx(gamma_valid, Pi_ls_valid * 100, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 5);
% hold on;
% xline(gamma_opt, 'r--', 'LineWidth', 1.5, ...
%     'Label', sprintf('\\gamma^*=%.2e', gamma_opt));
% yline(100, 'k--', 'LineWidth', 1, 'Label', '100% (LS baseline)');
% hold off;
% xlabel('$\gamma$',           'Interpreter','latex', 'FontSize', 12);
% ylabel('$\Pi(\gamma)$ (\%)', 'Interpreter','latex', 'FontSize', 12);
% title('$J_\gamma(b_\gamma) / J(b_{ls})$ (\%)', ...
%     'Interpreter','latex', 'FontSize', 12);
% grid on;
% 
% subplot(1,3,3);
% % Combined: n_active vs Pi — the Pareto frontier
% semilogx(n_active_valid, Pi_ls_valid * 100, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 5);
% hold on;
% scatter(n_active_valid(elbow_idx), Pi_ls_valid(elbow_idx)*100, ...
%     100, 'r', 'filled', 'DisplayName', sprintf('\\gamma^*: %d modes', ...
%     n_active_valid(elbow_idx)));
% hold off;
% xlabel('Active modes',       'FontSize', 12);
% ylabel('$\Pi(\gamma)$ (\%)', 'Interpreter','latex', 'FontSize', 12);
% title('Pareto: Modes vs Loss', 'FontSize', 12);
% legend('Location','best');
% grid on;
% 
% sgtitle(sprintf('Sparsity-Promoting Selection  |  \\gamma^* = %.2e  |  %d modes  |  \\Pi = %.2f%%', ...
%     gamma_opt, n_active_valid(elbow_idx), Pi_ls_valid(elbow_idx)*100), ...
%     'FontSize', 12);



%% ============================================================
% ADMM SOLVER
%% ============================================================
function b = sparsifyDMD(P, q, Pl, gamma, n)
% ADMM for L1-regularized amplitude:
% min (1/2)b'Pb - Re(q'b) + gamma*||b||_1
    rho   = 1;
    b     = zeros(n,1);
    beta  = zeros(n,1);
    lam   = zeros(n,1);
    maxIt = 1000;
    tol   = 1e-6;
    for iter = 1:maxIt
        b_old = b;
        b     = (Pl') \ ((Pl) \ (q + rho*(beta - lam)));
        beta  = max(0, abs(b+lam) - gamma/rho) .* sign(b+lam);
        lam   = lam + b - beta;
        if norm(b - b_old) < tol; break; end
    end
    b = beta;
end
