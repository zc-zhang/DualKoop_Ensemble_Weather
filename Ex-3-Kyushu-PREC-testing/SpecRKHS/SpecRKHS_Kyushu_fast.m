
%% Koopman/PF (dual RKHS ResDMD) analysis of ensemble Kyushu heavy data
%% Restructured one-to-one against the reference:
%% https://github.com/GustavConradie1/SpecRKHS/blob/main/examples/sealevels.m
rng(0)
T0 = tic;
tk = @(msg) fprintf('[%6.1fs] %s\n', toc(T0), msg);

%% =====================================================================
%  LOAD DATA
%  =====================================================================
% load('EnPREC_all_col_raw_72hrs.mat')
% Xa_all = EnVfull_Xa_36;
% Ya_all = EnVfull_Ya_36;

p         = size(Xa_all, 1);
M         = 100;
T_snap    = size(Xa_all, 2) / M;
grid_size = [140, 149];
nLAND     = (1:p).';

maxx  = 65;
steps = T_snap - maxx;

Xa_3D = reshape(Xa_all, p, T_snap, M);
Ya_3D = reshape(Ya_all, p, T_snap, M);

X = reshape(Xa_3D(:, 1:maxx, :), p, maxx * M);
Y = reshape(Ya_3D(:, 1:maxx, :), p, maxx * M);
n_train = maxx * M;

%% =====================================================================
%  GEO DATA + COASTLINE (cropped to the domain, merged into ONE line)
%  =====================================================================
lat_sub = weatherDat2021AUG_ensemble0.dat_lat;
lon_sub = weatherDat2021AUG_ensemble0.dat_lon;
lonlim  = [min(lon_sub(:)) max(lon_sub(:))];
latlim  = [min(lat_sub(:)) max(lat_sub(:))];

shapefile_path = 'D:\Susuki Lab\Testing_Code\data-weather\Data_250401\ne_10m_coastline\ne_10m_coastline.shp';
coast = shaperead(shapefile_path, 'BoundingBox', [lonlim(1) latlim(1); lonlim(2) latlim(2)]);
cx = []; cy = [];
for k = 1:numel(coast)
    cx = [cx, coast(k).X(:).', NaN]; %#ok<AGROW>
    cy = [cy, coast(k).Y(:).', NaN]; %#ok<AGROW>
end

% Kumamoto (32.803N, 130.708E)
targetLat = 32.803; targetLon = 130.708;
[~, Kumamoto_Locatidx] = min((lat_sub(:) - targetLat).^2 + (lon_sub(:) - targetLon).^2);
kuma_lon = lon_sub(Kumamoto_Locatidx);
kuma_lat = lat_sub(Kumamoto_Locatidx);
fprintf('Kumamoto index: %d\n', Kumamoto_Locatidx);
tk('geo data ready');

%% =====================================================================
%  KERNEL MATRICES
%  =====================================================================
ker = @(x, t) kernel_matern(x, t);
[G, A, R] = generate_matrices_kernelized(X, Y, ker);
disp('size(G):'); disp(size(G))
tk('G, A, R done');

%% =====================================================================
%  VERIFIED EIGENVALUES
%  =====================================================================
r_eig = 360;
num   = 36;
[Lambda_res, F_res, Lambda, F, res, res_verif, idx, W, W_res] = ...
    verified_eigenvalues(G, A, R, num, r_eig);
Lambda_res = Lambda_res(:);  Lambda = Lambda(:);
disp('length(idx):'); disp(length(idx))
tk('verified eigenvalues done');

%% =====================================================================
%  SHARED PRECOMPUTATION (done ONCE, not per member / per mode)
%  =====================================================================
X0  = squeeze(Xa_3D(:, maxx, :));                 % p x M  (start states)
K0  = kmat(X0, X);                                % M x n_train (replaces 2 scalar loops)
GF  = G * F_res;                                  % n_train x k
mean_state_train = mean(X, 1).';                  % n_train x 1

% one least-squares solve for ALL right-hand sides
sol       = GF \ [K0.', mean_state_train];
coefs_all = sol(:, 1:M).';                        % M x k
mode_res  = sol(:, M+1);                          % k x 1
tk('kernel vectors + LS solves done');

%% =====================================================================
%  PLOT: TOP 9 DUAL SPATIAL MODES
%  =====================================================================
mode_indices = [1,3,5,7,9,11,13,15,17];
mode_indices = mode_indices(mode_indices <= min([length(Lambda_res), size(F_res,2), length(res_verif)]));

% Phi for all modes at once:  ((G*F_j)\X.').'  ==  X*conj(g)/(g'*g)
GFm     = GF(:, mode_indices);
Phi_all = (X * conj(GFm)) ./ sum(abs(GFm).^2, 1);        % p x nModes

fig = figure('Color','w','Position',[100 100 1200 900]);
colormap(brighten(redblueTecplot(21), -0.55));
pp = 1;
for jj = 1:numel(mode_indices)
    j   = mode_indices(jj);
    L   = Lambda_res(j);
    r_j = res_verif(j);

    u = Phi_all(:, jj);
    u = real(u * exp(1i * mean(angle(u))));

    v = NaN(size(lon_sub));
    v(nLAND) = u(:);
    mode_map = abs(v);

    subplot(3, 3, pp);
    h = pcolor(lon_sub, lat_sub, mode_map);
    h.EdgeColor = 'none';
    shading interp; axis equal tight;
    set(gca, 'YDir', 'normal');

    c_max = max(mode_map(:), [], 'omitnan');
    c_min = min(mode_map(:), [], 'omitnan');
    if ~isnan(c_min) && ~isnan(c_max) && c_min ~= c_max
        clim([c_min, c_max]);
    end
    colorbar;

    hold on;
    plot(cx, cy, 'k-', 'LineWidth', 1.2);              % ONE call
    plot(kuma_lon, kuma_lat, 'p', 'MarkerSize', 12, ...
        'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
    hold off;

    title(sprintf('$V_{%d}$: $|\\lambda|=%.3f$, Res$=%.3f$', j, abs(L), r_j), ...
        'Interpreter', 'latex', 'FontSize', 10);
    xlabel('Lon ($^\circ$E)', 'Interpreter', 'latex');
    ylabel('Lat ($^\circ$N)', 'Interpreter', 'latex');
    xlim(lonlim); ylim(latlim);
    pp = pp + 1;
end
sgtitle('Dual Koopman spatial modes $V_j$', ...
    'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
exportgraphics(fig, 'top_9_dual_spatial_modes.pdf', ...
    'ContentType', 'image', 'Resolution', 300, 'BackgroundColor', 'none');
tk('mode figure done');

%% =====================================================================
%  PLOT SPURIOUS AND VERIFIED EIGENVALUES
%  =====================================================================
figure
scatter(angle(Lambda), log(abs(Lambda)), 200, res, '.', 'LineWidth', 1);
hold on
scatter(angle(Lambda_res), log(abs(Lambda_res)), 900, res_verif, '.', 'LineWidth', 1);
box on
clim([0, 0.04])
colormap(jet); colorbar
xlabel('$\mathrm{arg}(\lambda)$', 'interpreter', 'latex', 'fontsize', 18)
ylabel('$\mathrm{log}(|\lambda|)$', 'interpreter', 'latex', 'fontsize', 18)
title('Eigenvalues for ensemble data', 'interpreter', 'latex', 'fontsize', 18)
ax = gca; ax.FontSize = 18; axis([-pi pi -20*10^(-1) 10^(-1)])
for k = -5:1:5
    plot(k*pi/6*ones(22,1), -20*10^(-1):10^(-1):10^(-1), '--', 'Color', 'black')
end
xticks([-pi -5*pi/6 -4*pi/6 -3*pi/6 -2*pi/6 -pi/6 0 pi/6 2*pi/6 3*pi/6 4*pi/6 5*pi/6 pi])
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
xticklabels({'$-\pi$','','$-2\pi/3$','','$-\pi/3$','','$0$','','$\pi/3$','','$2\pi/3$','','$\pi$'})
xtickangle(30)
exportgraphics(gcf, 'ensemble_evals_angle.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

figure;
scatter(real(Lambda), imag(Lambda), 200, res, '.', 'LineWidth', 1);
hold on;
scatter(real(Lambda_res), imag(Lambda_res), 700, res_verif, '.', 'LineWidth', 1);
plot(cos(0:0.01:2*pi), sin(0:0.01:2*pi), '-k')
axis equal
axis([-1.15,1.15,-1.15,1.15])
clim([0,1])
colormap(jet); colorbar
xlabel('$\mathrm{Re}(\lambda)$', 'interpreter', 'latex', 'fontsize', 18)
ylabel('$\mathrm{Im}(\lambda)$', 'interpreter', 'latex', 'fontsize', 18)
title('Eigenvalues', 'interpreter', 'latex', 'fontsize', 18)
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
ax = gca; ax.FontSize = 18; box on;
exportgraphics(gcf, 'KE_evals.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')
tk('eigenvalue figures done');

%% =====================================================================
%  PREDICTIONS: SpecRKHS-Obs (KMD) and kEDMD  -- no per-member big solves
%  =====================================================================
x_kmd_all   = zeros(p, steps, M);
x_kedmd_all = zeros(p, steps, M);

% ---- KMD ----
FX     = F_res.' * X.';                            % k x p, once
lamPow = (Lambda_res .^ (1:steps)).';              % steps x k
for m = 1:M
    x_kmd_all(:, :, m) = real((coefs_all(m, :) .* lamPow) * FX).';
end
tk('KMD forecasts done');

% ---- kEDMD ----
% Original: mode_full = ([G;g0]*W) \ [X, x0].'  (LS with p right-hand sides,
% repeated for every member).  Equivalent cheap form: economy QR of G*W once,
% then only a small ((k+1) x k) QR per member.  Same LS solution.
GW      = G * W;
[Qw,Rw] = qr(GW, 0);
XQ      = X * Qw;                                  % p x k
QtMean  = Qw' * mean_state_train;                  % k x 1
PSI0    = K0 * W;                                  % M x k
kW      = size(W, 2);

mean_obs_kedmd_all = zeros(steps+1, M);
for m = 1:M
    x0_m = X0(:, m);
    psi0 = PSI0(m, :);
    [Q2, R2] = qr([Rw; psi0], 0);

    Dc = psi0' .* (Lambda .^ (1:steps));           % k x steps
    Z  = Q2 * (R2' \ Dc);
    x_kedmd_all(:, :, m) = real(XQ * Z(1:kW, :) + x0_m * Z(kW+1, :));

    % domain-mean observable via kEDMD
    mo = R2 \ (Q2' * [QtMean; mean(x0_m)]);        % k x 1
    Dobs = psi0.' .* (conj(Lambda) .^ (0:steps));  % k x (steps+1)
    mean_obs_kedmd_all(:, m) = real(Dobs.' * mo);
end
tk('kEDMD forecasts done');

%% =====================================================================
%  DMD (pooled training) -- forecast done in reduced coordinates
%  =====================================================================
[U, Ssv, ~] = svd(X, 'econ');
sv = diag(Ssv);
r_dmd = nnz(sv > max(size(X)) * eps(sv(1)));       % same as rank(), no 2nd SVD
r_dmd_max = Inf;   % <-- set e.g. 1000 for a big speed-up if sv decays fast
r_dmd = min(r_dmd, r_dmd_max);
U = U(:, 1:r_dmd);
tk('SVD done');

PXs = X.' * U;
PYs = Y.' * U;
Kd  = PXs \ PYs;
[Wdmd, LAM, ~] = eig(Kd, 'vector');
tk('DMD eig done');

Cdmd = Wdmd \ (U.' * X0);                          % r x M   (== Phi_dmd \ x0_m)
x_dmd_all = zeros(p, steps, M);
for t = 1:steps
    x_dmd_all(:, t, :) = reshape(real(U * (Wdmd * (Cdmd .* (LAM .^ t)))), p, 1, M);
end
tk('DMD forecasts done');

%% =====================================================================
%  RELATIVE FORECAST ERRORS
%  =====================================================================
er1_all = zeros(M, steps);   % SpecRKHS-Obs
er2_all = zeros(M, steps);   % DMD
er3_all = zeros(M, steps);   % kEDMD
for m = 1:M
    real_data_m = Ya_3D(:, maxx + (1:steps), m);
    den = sum(abs(real_data_m).^2, 1);
    er1_all(m, :) = sum(abs(x_kmd_all(:, :, m)   - real_data_m).^2, 1) ./ den;
    er2_all(m, :) = sum(abs(x_dmd_all(:, :, m)   - real_data_m).^2, 1) ./ den;
    er3_all(m, :) = sum(abs(x_kedmd_all(:, :, m) - real_data_m).^2, 1) ./ den;
end
tk('errors computed');

N_steps = size(er1_all, 2);
t_vec   = 1:N_steps;
x_poly  = [t_vec, fliplr(t_vec)];

%% ---- Figure: all members + means (linear) ----
figure
hold on
plot_members(er2_all, [0.4 0.4 0.4], 0.2);
plot_members(er3_all, [0.9 0.4 0.2], 0.2);
plot_members(er1_all, [0.2 0.4 0.8], 0.2);
plot(mean(er2_all, 1), 'k--', 'linewidth', 2, 'DisplayName', 'DMD (mean)')
plot(mean(er3_all, 1), 'Color', [0.8500 0.3250 0.0980], 'linewidth', 2, 'DisplayName', 'KeDMD (mean)')
plot(mean(er1_all, 1), 'b', 'linewidth', 2, 'DisplayName', 'DualKoop-RKHS (mean)')
grid on; box on;
title('Relative forecast errors comparison', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend('interpreter', 'latex', 'fontsize', 16, 'location', 'best')
exportgraphics(gcf, 'ensemble_weather_error.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% ---- Figure: all members + means (log) with zoom inset ----
c_dmd   = [0.4, 0.4, 0.4];
c_kedmd = [0.8500, 0.3250, 0.0980];
c_dual  = [0.2, 0.4, 0.8];

figure('Color', 'w', 'Position', [100, 100, 900, 650]);
hold on;
plot_members(er2_all, c_dmd,   0.15);
plot_members(er3_all, c_kedmd, 0.15);
plot_members(er1_all, c_dual,  0.15);
plot(mean(er2_all, 1), '--', 'Color', c_dmd,   'LineWidth', 2, 'DisplayName', 'DMD (mean)')
plot(mean(er3_all, 1), '-',  'Color', c_kedmd, 'LineWidth', 2, 'DisplayName', 'KeDMD (mean)')
plot(mean(er1_all, 1), '-.', 'Color', c_dual,  'LineWidth', 2, 'DisplayName', 'DualKoop-RKHS (mean)')
set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 14);
grid on; box on;
title('Relative Forecast Errors Comparison', 'Interpreter', 'latex', 'FontSize', 18);
xlabel('Lead Time (Snapshots)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Relative Forecast Error', 'Interpreter', 'latex', 'FontSize', 18);
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northwest');

axes('Position', [0.48, 0.22, 0.40, 0.38]);
hold on;
plot_members(er3_all, c_kedmd, 0.2);
plot_members(er1_all, c_dual,  0.2);
plot(mean(er3_all, 1), '-',  'Color', c_kedmd, 'LineWidth', 2);
plot(mean(er1_all, 1), '-.', 'Color', c_dual,  'LineWidth', 2);
set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 11);
grid on; box on;
title('\textbf{Zoomed: KeDMD vs DualKoop}', 'Interpreter', 'latex', 'FontSize', 11);
xlim([1, size(er1_all, 2)]);
ylim([min([er1_all(:); er3_all(:)]) * 0.9, max([er1_all(:); er3_all(:)]) * 1.1]);
exportgraphics(gcf, 'ensemble_weather_error_1.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%% ---- Figure: shaded min-max range ----
min1 = min(er1_all, [], 1); max1 = max(er1_all, [], 1);
min2 = min(er2_all, [], 1); max2 = max(er2_all, [], 1);
min3 = min(er3_all, [], 1); max3 = max(er3_all, [], 1);

figure
hold on
fill(x_poly, [max2, fliplr(min2)], [0.0000 0.4470 0.7410], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off')
fill(x_poly, [max3, fliplr(min3)], [0.8500 0.3250 0.0980], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off')
fill(x_poly, [max1, fliplr(min1)], [0.9290 0.6940 0.1250], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off')
plot(mean(er2_all, 1), '--', 'Color', [0.0000 0.4470 0.7410], 'LineWidth', 2, 'DisplayName', 'DMD (mean)')
plot(mean(er3_all, 1), '-',  'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, 'DisplayName', 'KeDMD (mean)')
plot(mean(er1_all, 1), '-.', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2, 'DisplayName', 'DualKoop-RKHS (mean)')
set(gca, 'YScale', 'log');
grid on; box on;
title('Relative forecast errors comparison', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend('interpreter', 'latex', 'fontsize', 16, 'location', 'best')
exportgraphics(gcf, 'ensemble_shallow_weather_error.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% ---- Figure: shaded range + zoom inset ----
c_dmd   = [0.0000, 0.4470, 0.7410];
c_kedmd = [0.8500, 0.3250, 0.0980];
c_dual  = [0.9290, 0.6940, 0.1250];

figure('Color', 'w', 'Position', [100, 100, 900, 650]);
hold on;
fill(x_poly, [max2, fliplr(min2)], c_dmd,   'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(x_poly, [max3, fliplr(min3)], c_kedmd, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(x_poly, [max1, fliplr(min1)], c_dual,  'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(mean(er2_all, 1), '--', 'Color', c_dmd,   'LineWidth', 2, 'DisplayName', 'DMD (mean)');
plot(mean(er3_all, 1), '-',  'Color', c_kedmd, 'LineWidth', 2, 'DisplayName', 'KeDMD (mean)');
plot(mean(er1_all, 1), '-.', 'Color', c_dual,  'LineWidth', 2, 'DisplayName', 'DualKoop-RKHS (mean)');
set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 14);
grid on; box on;
title('Relative Forecast Errors Comparison', 'Interpreter', 'latex', 'FontSize', 18);
xlabel('Lead Time (Snapshots)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Relative Forecast Error', 'Interpreter', 'latex', 'FontSize', 18);
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northwest');

axes('Position', [0.50, 0.22, 0.38, 0.38]);
hold on;
fill(x_poly, [max3, fliplr(min3)], c_kedmd, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(x_poly, [max1, fliplr(min1)], c_dual,  'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(mean(er3_all, 1), '-',  'Color', c_kedmd, 'LineWidth', 2);
plot(mean(er1_all, 1), '-.', 'Color', c_dual,  'LineWidth', 2);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
grid on; box on;
title('\textbf{Zoomed: KeDMD vs DualKoop}', 'Interpreter', 'latex', 'FontSize', 11);
xlim([1, N_steps]);
ylim([min([min1, min3]) * 0.9, max([max1, max3]) * 1.1]);
exportgraphics(gcf, 'ensemble_shallow_weather_error_1.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%% ---- Bar plots: overall RMSE ----
rmse_values = [sqrt(mean(er2_all(:).^2)); sqrt(mean(er3_all(:).^2)); sqrt(mean(er1_all(:).^2))];
c_dmd   = [0.0000 0.4470 0.7410];
c_kedmd = [0.8500 0.3250 0.0980];
c_spec  = [0.9290 0.6940 0.1250];

figure('Color', 'w');
b = bar(1:3, rmse_values, 'FaceColor', 'flat', 'BarWidth', 0.5);
b.CData(1,:) = c_dmd; b.CData(2,:) = c_kedmd; b.CData(3,:) = c_spec;
set(gca, 'YScale', 'log');
set(gca, 'XTick', 1:3, 'XTickLabel', {'DMD', 'KeDMD', 'DualKoop-RKHS'});
grid on; box on;
ylabel('Overall RMSE', 'interpreter', 'latex', 'fontsize', 18);
title('Overall Forecast Error Comparison (RMSE)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'FontSize', 14);
exportgraphics(gcf, 'ensemble_overall_rmse_bar.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

figure('Color', 'w', 'Position', [100, 100, 800, 600]);
b = bar(1:3, rmse_values, 'FaceColor', 'flat', 'BarWidth', 0.5);
b.CData(1,:) = c_dmd; b.CData(2,:) = c_kedmd; b.CData(3,:) = c_spec;
set(gca, 'TickLabelInterpreter', 'latex', 'XTick', 1:3, ...
    'XTickLabel', {'\textrm{DMD}', '\textrm{KeDMD}', '\textrm{DualKoop}'});
grid on; box on;
ylabel('Overall RMSE', 'Interpreter', 'latex', 'FontSize', 18);
title('Overall Forecast Error Comparison (RMSE)', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 14);

axes('Position', [0.50, 0.32, 0.38, 0.38]);
b_inset = bar(2:3, rmse_values(2:3), 'FaceColor', 'flat', 'BarWidth', 0.5);
b_inset.CData(1,:) = c_kedmd; b_inset.CData(2,:) = c_spec;
set(gca, 'TickLabelInterpreter', 'latex', 'XTick', 2:3, ...
    'XTickLabel', {'\textrm{KeDMD}', '\textrm{DualKoop-RKHS}'});
title('\textbf{Zoomed (KeDMD vs DualKoop)}', 'Interpreter', 'latex', 'FontSize', 11);
grid on; box on;
set(gca, 'FontSize', 11);
exportgraphics(gcf, 'ensemble_overall_rmse_bar_1.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
tk('error figures done');

%% =====================================================================
%  DOMAIN-MEAN OBSERVABLE (all members vectorized)
%  =====================================================================
pow = (Lambda_res .^ (0:steps)).';                 % (steps+1) x k

% KMD prediction of the observable
v_obs = F_res.' * (GF * mode_res);                 % == F_res.'*G*F_res*mode_res
mean_obs_kmd_all = real(pow * (coefs_all.' .* v_obs));           % (steps+1) x M

% exact domain mean over the forecast window
mean_obs_exact_all = squeeze(mean(Ya_3D(:, maxx-1 + (1:steps+1), :), 1));   % (steps+1) x M

% exact kernel-space observable: all future points x all members in ONE kmat call
Yfut = reshape(Ya_3D(:, maxx-1 + (1:steps+1), :), p, (steps+1) * M);
Kf   = kmat(Yfut, X);                              % ((steps+1)M) x n_train
obs_exact_all = reshape(Kf * (F_res * mode_res), steps+1, M);

% K_{x0} prediction error
FK          = F_res.' * K0.';                                    % k x M
Kx0_predict = real(pow * (coefs_all.' .* FK));                   % (steps+1) x M
Kx0_future  = zeros(steps+1, M);
for m = 1:M
    cols = (m-1)*(steps+1) + (1:steps+1);
    Kx0_future(:, m) = kmat(X0(:, m), Yfut(:, cols)).';
end
Kx0_er_all = abs(Kx0_future - Kx0_predict).^2 ./ abs(Kx0_future).^2;
tk('domain-mean observables done');

%% ---- Mean forecast comparison ----
figure
hold on
p1 = plot(0:steps, log(mean(mean_obs_kedmd_all, 2)), 'linewidth', 3, 'color', [0.8500 0.3250 0.0980]);
p2 = plot(0:steps, log(mean(mean_obs_kmd_all, 2)), '-.', 'linewidth', 3, 'color', [0.9290 0.6940 0.1250]);
p3 = plot(0:steps, log(mean(mean_obs_exact_all, 2)), '--', 'linewidth', 3, 'color', [0 0.4470 0.7410]);
grid on; box on;
title('Mean forecast comparison', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Domain-Mean Prediction', 'interpreter', 'latex', 'fontsize', 18)
legend([p3 p1 p2], {'Exact', 'KeDMD', 'DualKoop-RKHS'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'best')
exportgraphics(gcf, 'ensemble_prediction_mean.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% ---- Domain-mean relative error ----
figure
er_obs = mean(abs(obs_exact_all - mean_obs_kmd_all).^2 ./ abs(obs_exact_all).^2, 2);
semilogy(0:steps, er_obs, 'linewidth', 3, 'color', [0.9290 0.6940 0.1250])
grid on
title('Domain-mean relative forecast error', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend({'DualKoop-RKHS'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'best')
exportgraphics(gcf, 'ensemble_error_mean.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% ---- Kernel-function relative error ----
figure
semilogy(0:steps, mean(Kx0_er_all, 2), 'linewidth', 3, 'color', [0.9290 0.6940 0.1250])
grid on
title('Relative forecast error for kernel function', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend({'DualKoop-RKHS'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'best')
exportgraphics(gcf, 'ensemble_prediction_kernel_error.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')
tk('ALL DONE');

%% =====================================================================
%  LOCAL FUNCTIONS
%  =====================================================================
function ker = kernel_matern(x, t)
    % original scalar/column version (used by generate_matrices_kernelized)
    sigma = 1/20000;
    r = vecnorm(x - t);
    ker = zeros(1, size(x, 2));
    ker(r>0)  = (sigma*r(r>0)).^(3/2) .* besselk(-3/2, sigma*r(r>0));
    ker(r==0) = sqrt(pi/2);
end

function K = kmat(A, B)
    % Vectorized Matern kernel between all columns of A (p x na) and B (p x nb)
    % -> na x nb.  Same formula as kernel_matern.
    sigma = 1/20000;
    D2 = sum(A.^2, 1).' + sum(B.^2, 1) - 2 * (A.' * B);
    r  = sqrt(max(D2, 0));
    K  = sqrt(pi/2) * ones(size(r));
    nz = r > 0;
    z  = sigma * r(nz);
    K(nz) = z.^(3/2) .* besselk(-3/2, z);
end

function plot_members(E, col, a)
    % draw every row of E (members x steps) as a translucent line, hidden from legend
    h = plot(E.', 'Color', [col, a], 'LineWidth', 1);
    set(h, 'HandleVisibility', 'off');
end