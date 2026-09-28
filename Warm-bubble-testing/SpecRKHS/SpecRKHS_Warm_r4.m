%% Koopman/PF (dual RKHS ResDMD) analysis of ensemble weather-simulation data
%% Restructured one-to-one against the reference:
%% https://github.com/GustavConradie1/SpecRKHS/blob/main/examples/sealevels.m
rng(0)

%% =====================================================================
%  LOAD DATA
%  =====================================================================
% Expect EnVfull_Xa_36, EnVfull_Ya_36 in the workspace (or load them here,
% mirroring sealevels.m's load('DATA_north.mat') etc.)
% load('EnVfull_data.mat')

Xa_all = EnVfull_Xa_36;
Ya_all = EnVfull_Ya_36;

p         = size(Xa_all, 1);   % state dimension
M         = 10;                % number of ensemble members
T_snap    = size(Xa_all, 2) / M;  % time snapshots per member
grid_size = [97, 40];          % [rows(z), cols(y)] for spatial plots
nLAND     = (1:p).';           % identity mask (no land/ocean exclusion here)

%% =====================================================================
%  SEPARATE INTO TRAINING AND TEST DATA SETS
%  (sealevels.m: x = DATA(:,min:maxx-1); y = DATA(:,min+1:maxx))
%  Here the same one-step-delay pairing is done *within each member* and
%  the pairs are pooled across members.
%  =====================================================================
maxx  = 30;             % training snapshots used per member
steps = T_snap - maxx;  % forecast horizon (remaining snapshots)

Xa_3D = reshape(Xa_all, p, T_snap, M);
Ya_3D = reshape(Ya_all, p, T_snap, M);

X = reshape(Xa_3D(:, 1:maxx, :), p, maxx * M);   % pooled "x"
Y = reshape(Ya_3D(:, 1:maxx, :), p, maxx * M);   % pooled "y"
n_train = maxx * M;

%% =====================================================================
%  COMPUTE MATRICES
%  =====================================================================
ker = @(x, t) kernel_matern(x, t);
% 1. Parse inputs or specify kernel type
%kernel_type = "Gaussian"; % or p.Results.type
% 
% % 2. Get the handle
% ker = get_kernel(kernel_type, X);

[G, A, R] = generate_matrices_kernelized(X, Y, ker);
disp('size(G):'); disp(size(G))

%% =====================================================================
%  COMPUTE VERIFIED EIGENVALUES
%  sealevels.m keeps 5 eigenpairs for a scalar sea-level field; the
%  ensemble state here is higher-dimensional, so keep more.
%  =====================================================================
r =280;
num=17;
[Lambda_res, F_res, Lambda, F, res, res_verif, idx, W, W_res] = ...
    verified_eigenvalues(G, A, R, num,r);
disp('length(idx):'); disp(length(idx))

%% =====================================================================
%  COMPUTE PERRON-FROBENIUS / KOOPMAN MODES AND PLOT (single-mode panels)
%  Mirrors sealevels.m's `for idx=[2 4]` loop: pick a couple of
%  representative verified modes and plot each on its own figure.
%  =====================================================================
n_modes_available = length(Lambda_res);
%single_mode_idx = unique(round(linspace(1, n_modes_available, min(17, n_modes_available))));
single_mode_idx =[1,3,5,7];
for j = single_mode_idx
    L = Lambda_res(j);
    F_j = F_res(:, j);
    r_j = res_verif(j);

    Phi = ((G * F_j) \ (X.')).';

    u = Phi;
    u = real(u * exp(1i * mean(angle(u))));

    v = zeros(grid_size(1) * grid_size(2), 1) + NaN;
    v(nLAND) = u(:);
    v = reshape(v, grid_size);

    figure
    imagesc(data.y, data.z, abs(v), 'AlphaData', ~isnan(v))
    colormap(brighten(redblueTecplot(21), -0.55));
    colorbar
    %clim([mean(u(:)) - 2*std(u(:)), mean(u(:)) + 2*std(u(:))])
    set(gca, 'Color', [1,1,1]*0.6,'FontSize',16)
    axis xy; axis equal
    xlim([0 2e4]); ylim([0 2e4])
   % grid on; 
    box on
    title(sprintf('Res$=$%f, $\\lambda=%.2f%+.2fi$', abs(r_j), real(L), imag(L)), ...
        'interpreter', 'latex', 'fontsize', 17)
    exportgraphics(gcf, sprintf('ensemble_mode_%d.pdf', j), ...
        'ContentType', 'vector', 'BackgroundColor', 'none')
end


%% Plot 3 x 3 outplots of spatial figures 

% =====================================================================
n_modes_available = length(Lambda_res);

% Select 9 evenly spaced indices up to 17 (or max available) for a 3x3 grid
max_modes = min(17, n_modes_available);
selected_indices = unique(round(linspace(1, max_modes, 9)));

% Create a single figure window for the 3x3 grid
figure('Position', [100, 100, 1100, 950]); 

% Use 'tight' spacing to reduce gaps between columns and rows
tiledlayout(3, 3, 'TileSpacing', 'tight', 'Padding', 'tight');

for k = 1:length(selected_indices)
    j = selected_indices(k);
    
    L = Lambda_res(j);
    F_j = F_res(:, j);
    r_j = res_verif(j);
    Phi = ((G * F_j) \ (X.')).';
    u = Phi;
    u = real(u * exp(1i * mean(angle(u))));
    
    v = zeros(grid_size(1) * grid_size(2), 1) + NaN;
    v(nLAND) = u(:);
    v = reshape(v, grid_size);
    
    % Place in the next tile of the 3x3 grid
    nexttile;
    
    imagesc(data.y, data.z, abs(v), 'AlphaData', ~isnan(v))
    colormap(gca, brighten(redblueTecplot(21), -0.55));
    colorbar
    
    set(gca, 'Color', [1,1,1]*0.6,'FontSize',11)
    axis xy; axis equal
    xlim([0 2e4]); ylim([0 2e4])
    box on
    
    % Include \lambda_{j} with the explicit index in the title
    title(sprintf('Res$=$%.4f, $\\lambda_{%d}=%.2f%+.2fi$', abs(r_j), j, real(L), imag(L)), ...
        'interpreter', 'latex', 'fontsize', 8.5)
end

% Export the entire 3x3 figure as a single vector PDF
exportgraphics(gcf, 'ensemble_modes_3x3_grid.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none')


%% =====================================================================
%  PLOT SPURIOUS AND VERIFIED EIGENVALUES
%  =====================================================================
figure;
scatter(angle(Lambda), log(abs(Lambda)), 150, res, '.', 'LineWidth', 1);
hold on
scatter(angle(Lambda_res), log(abs(Lambda_res)), 800, res_verif, '.', 'LineWidth', 1);
box on
clim([0, 0.04])
colormap('turbo'); colorbar
xlabel('$\mathrm{arg}(\lambda)$', 'interpreter', 'latex', 'fontsize', 18)
ylabel('$\mathrm{log}(|\lambda|)$', 'interpreter', 'latex', 'fontsize', 18)
title(['Eigenvalues for ensemble data'], 'interpreter', 'latex', 'fontsize', 18)
ax = gca; ax.FontSize = 18; 

% Option 1: Let MATLAB auto-scale the y-axis automatically
xlim([-pi pi]);
 ylim([-0.5, 0.5]); % <-- Option 2: Alternatively, set a much wider manual limit if needed

for k = -5:1:5
    % Adjusted the vertical lines to cover a broader y-range if you use manual limits
    plot(k*pi/6*ones(200,1), linspace(ax.YLim(1), ax.YLim(2), 200), '--', 'Color', [0.7 0.7 0.7])
end

xticks([-pi -5*pi/6 -4*pi/6 -3*pi/6 -2*pi/6 -pi/6 0 pi/6 2*pi/6 3*pi/6 4*pi/6 5*pi/6 pi])
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
xticklabels({'$-\pi$','','$-2\pi/3$','','$-\pi/3$','','$0$','','$\pi/3$','','$2\pi/3$','','$\pi$'})
xtickangle(30)
exportgraphics(gcf, 'ensemble_evals_angle.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% ==================================Eigvalues 
figure;
scatter(real(Lambda),imag(Lambda),300,res,'.','LineWidth',1);
hold on
scatter(real(Lambda_res),imag(Lambda_res),600,res_verif,'.','LineWidth',1);
hold on
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
axis equal
axis([-1.15,1.15,-1.15,1.15])
clim([0,1])
load('cmap.mat')
colormap(jet); colorbar
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
title(sprintf('Eigenvalues'),'interpreter','latex','fontsize',18)
%ax=gca; ax.FontSize=18; box on;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
ax = gca; 
ax.FontSize = 18; box on;
% Vector Export (from Part 2)
exportgraphics(gcf, 'KE_evals.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%====================

%% =====================================================================
%  PSEUDOSPECTRA  (same as sealevels.m; expensive -- optional)
%  =====================================================================
% pts   = 100;
% x_pts = linspace(-1.2, 1.2, pts);
% y_pts = linspace(-0.02, 1.2, pts/2);
% z_pts = kron(x_pts, ones(length(y_pts),1)) + 1i*kron(ones(1,length(x_pts)), y_pts(:));
% z_pts = z_pts(:);
% res_pspec = pseudospectra(G, A, R, z_pts);
% res_pspec_rs = reshape(res_pspec, length(y_pts), length(x_pts));
% figure; hold on; box on
% v = (10.^(-10:0.3:0));
% contourf(reshape(real(z_pts),length(y_pts),length(x_pts)), reshape(imag(z_pts),length(y_pts),length(x_pts)), log10(real(res_pspec_rs)), log10(v));
% contourf(reshape(real(z_pts),length(y_pts),length(x_pts)), -reshape(imag(z_pts),length(y_pts),length(x_pts)), log10(real(res_pspec_rs)), log10(v));
% cbh = colorbar; cbh.Ticks = log10(10.^(-2:1:0)); cbh.TickLabels = 10.^(-2:1:0);
% clim([-2,0]); reset(gcf); set(gca, 'YDir', 'normal'); colormap('parula'); axis equal;
% plot(sin(0:0.01:2*pi), cos(0:0.01:2*pi), '--', 'color', 'white');
% grid on
% title('Pseudospectrum of ensemble data', 'interpreter', 'latex', 'fontsize', 18)
% xlabel('$\mathrm{Re}(z)$', 'interpreter', 'latex', 'fontsize', 18)
% ylabel('$\mathrm{Im}(z)$', 'interpreter', 'latex', 'fontsize', 18)
% ax = gca; ax.FontSize = 18; axis([x_pts(1), x_pts(end), -y_pts(end), y_pts(end)])
% exportgraphics(gcf, 'ensemble_pspec.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% =====================================================================
%  COMPUTE PREDICTIONS USING SpecRKHS-Obs  (per member, generalizing
%  sealevels.m's single x0=DATA(:,maxx))
%  =====================================================================
x_kmd_all   = zeros(p, steps, M);
x_kedmd_all = zeros(p, steps, M);

for m = 1:M
    x0_m = Xa_3D(:, maxx, m);

    Kx0_vals = zeros(n_train, 1);
    for i = 1:n_train
        Kx0_vals(i) = ker(x0_m, X(:, i));
    end
    coefs_res = ((G * F_res) \ Kx0_vals).';
    x_kmd_all(:, :, m) = real((coefs_res .* (Lambda_res .^ (1:steps)).') * (F_res.' * X.')).';

    %% kEDMD predictions instead using KMD
    G_start = zeros(1, n_train);
    for i = 1:n_train
        G_start(i) = ker(x0_m, X(:, i));
    end
    mode_full = (([G; G_start] * W) \ ([X, x0_m].')).';
    psi0_full = G_start * W;
    x_kedmd_all(:, :, m) = real(transpose(transpose(psi0_full) .* (conj(Lambda).^(1:steps))) * mode_full.')';
end

%% =====================================================================
%  COMPARE TO DMD  (pooled training via SVD, per-member forecast from x0_m)
%  Note: sealevels.m instead builds a single companion-form predictor
%  from one long trajectory (PXr, PYr, c). With M separate members there
%  is no single "next snapshot" to append, so each member's forecast is
%  built from its own x0_m projected onto the pooled DMD modes -- the
%  natural per-member generalization of that approach.
%  =====================================================================
[U, S, ~] = svd(X, 'econ');
r = rank(S);
U = U(:, 1:r);

PXs = X' * U;
PYs = Y' * U;
K = PXs \ PYs;
[Wdmd, LAM, ~] = eig(K, 'vector');

Phi_dmd = U * Wdmd;   % p x r

x_dmd_all = zeros(p, steps, M);
for m = 1:M
    x0_m = Xa_3D(:, maxx, m);
    c_m = Phi_dmd \ x0_m;
    Time_Dynamics = c_m .* (LAM .^ (1:steps));
    x_dmd_all(:, :, m) = real(Phi_dmd * Time_Dynamics);
end

%% =====================================================================
%  PLOT RELATIVE FORECAST ERRORS  (per member, then mean +/- spread,
%  generalizing sealevels.m's single er1/er2/er3 curves)
%  =====================================================================
er1_all = zeros(M, steps);   % SpecRKHS-Obs
er2_all = zeros(M, steps);   % DMD
er3_all = zeros(M, steps);   % kEDMD

for m = 1:M
    real_data_m = Ya_3D(:, maxx + (1:steps), m);
    er1_all(m, :) = sum(abs(x_kmd_all(:, :, m)   - real_data_m).^2, 1) ./ sum(abs(real_data_m).^2, 1);
    er2_all(m, :) = sum(abs(x_dmd_all(:, :, m)   - real_data_m).^2, 1) ./ sum(abs(real_data_m).^2, 1);
    er3_all(m, :) = sum(abs(x_kedmd_all(:, :, m) - real_data_m).^2, 1) ./ sum(abs(real_data_m).^2, 1);
end

figure
hold on
for m = 1:M
    plot(er2_all(m,:), 'Color', [0.4 0.4 0.4 0.2], 'LineWidth', 1, 'HandleVisibility', 'off')
    plot(er3_all(m,:), 'Color', [0.9 0.4 0.2 0.2], 'LineWidth', 1, 'HandleVisibility', 'off')
    plot(er1_all(m,:), 'Color', [0.2 0.4 0.8 0.2], 'LineWidth', 1, 'HandleVisibility', 'off')
end
plot(mean(er2_all, 1), 'k--', 'linewidth', 2, 'DisplayName', 'DMD (mean)')
plot(mean(er3_all, 1), 'Color', [0.8500 0.3250 0.0980], 'linewidth', 2, 'DisplayName', 'KeDMD (mean)')
plot(mean(er1_all, 1), 'b', 'linewidth', 2, 'DisplayName', 'DualKoop-RKHS (mean)')
grid on; box on;
title('Relative forecast errors comparison', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend('interpreter', 'latex', 'fontsize', 16, 'location', 'best')
exportgraphics(gcf, 'ensemble_weather_error.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')


%% ------- Lots os lines of ensembles 
figure('Color', 'w', 'Position', [100, 100, 900, 650]);
hold on;
% Define colors
c_dmd = [0.4, 0.4, 0.4];
c_ked = [0.8500, 0.3250, 0.0980];
%c_dual = [0.2, 0.4, 0.8]; 
c_dual = [0.9290, 0.6940, 0.1250];
% --- 1. Main Plot (All Ensembles + Means) ---
M = size(er1_all, 1);
for m = 1:M
plot(er2_all(m,:), 'Color', [c_dmd, 0.15], 'LineWidth', 1, 'HandleVisibility', 'off')
plot(er3_all(m,:), 'Color', [c_ked, 0.15], 'LineWidth', 1, 'HandleVisibility', 'off')
plot(er1_all(m,:), 'Color', [c_dual, 0.15], 'LineWidth', 1, 'HandleVisibility', 'off')
end
plot(mean(er2_all, 1), '--', 'Color', c_dmd, 'LineWidth', 2, 'DisplayName', 'DMD (mean)')
plot(mean(er3_all, 1), '-', 'Color', c_ked, 'LineWidth', 2, 'DisplayName', 'KeDMD (mean)')
plot(mean(er1_all, 1), '-.', 'Color', c_dual, 'LineWidth', 2, 'DisplayName', 'DualKoop-RKHS (mean)')
% Main Axes Formatting
set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 25);
grid on; box on;
title('Relative Forecast Errors Comparison', 'Interpreter', 'latex', 'FontSize', 20);
xlabel('Lead Time (Snapshots)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Relative Forecast Error', 'Interpreter', 'latex', 'FontSize', 20);
legend('Interpreter', 'latex', 'FontSize', 16, 'Location', 'northwest');
% % --- 2. Inset Zoom Plot (KeDMD vs DualKoop-RKHS in Log Scale) ---
% % Position: [left, bottom, width, height]
% axes('Position', [0.48, 0.22, 0.40, 0.38]);
% hold on;
% % Plot ensemble lines for KeDMD and DualKoop only
% for m = 1:M
% plot(er3_all(m,:), 'Color', [c_ked, 0.2], 'LineWidth', 0.8, 'HandleVisibility', 'off');
% plot(er1_all(m,:), 'Color', [c_dual, 0.2], 'LineWidth', 0.8, 'HandleVisibility', 'off');
% end
% % Plot mean lines
% plot(mean(er3_all, 1), '-', 'Color', c_ked, 'LineWidth', 2);
% plot(mean(er1_all, 1), '-.', 'Color', c_dual, 'LineWidth', 2);
% % Inset Formatting (Log Scale)
% set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 11);
% grid on; box on;
% title('\textbf{Zoomed: KeDMD vs DualKoop}', 'Interpreter', 'latex', 'FontSize', 11);
% % Auto-tight y-limits around KeDMD and DualKoop data
% y_min_zoom = min([er1_all(:); er3_all(:)]);
% y_max_zoom = max([er1_all(:); er3_all(:)]);
% xlim([1, size(er1_all, 2)]);
% ylim([y_min_zoom * 0.9, y_max_zoom * 1.1]);
% --- 3. Export Vector Graphics PDF ---
exportgraphics(gcf, 'ensemble_weather_error_lines.pdf', ...
'ContentType', 'vector', 'BackgroundColor', 'none');




%%  ==============shallow plot the relative error
figure
hold on

% 1. Time axis length based on columns of er1_all
N_steps = size(er1_all, 2);
t_vec   = 1:N_steps;
x_poly  = [t_vec, fliplr(t_vec)];

% 2. Compute min-max shaded ranges across ensembles (dim 1)
min1 = min(er1_all, [], 1); max1 = max(er1_all, [], 1);
min2 = min(er2_all, [], 1); max2 = max(er2_all, [], 1);
min3 = min(er3_all, [], 1); max3 = max(er3_all, [], 1);

% 3. Plot shaded ranges (fill)
fill(x_poly, [max2, fliplr(min2)], [0.0000 0.4470 0.7410], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off')
fill(x_poly, [max3, fliplr(min3)], [0.8500 0.3250 0.0980], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off')
fill(x_poly, [max1, fliplr(min1)], [0.9290 0.6940 0.1250], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off')

% 4. Plot mean lines
plot(mean(er2_all, 1), '--', 'Color', [0.0000 0.4470 0.7410], 'LineWidth', 3, 'DisplayName', 'DMD (mean)')
plot(mean(er3_all, 1), '', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 3, 'DisplayName', 'KeDMD (mean)')
plot(mean(er1_all, 1), '-.', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 3, 'DisplayName', 'DualKoop-RKHS (mean)')

% 5. Log-scale y-axis formatting (replaces invalid log(...) calls)
set(gca, 'YScale', 'log','FontSize', 25);

grid on; box on;
title('Relative forecast errors comparison', 'fontsize', 20, 'interpreter', 'latex')
xlabel('Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend('interpreter', 'latex', 'fontsize', 16, 'location', 'best')
exportgraphics(gcf, 'ensemble_shallow_weather_error.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')


%% Additional Zoon plot of shallow plot
figure('Color', 'w', 'Position', [100, 100, 900, 650]);
hold on;
% --- 1. Data Preparation ---
N_steps = size(er1_all, 2);
t_vec = 1:N_steps;
x_poly = [t_vec, fliplr(t_vec)];
% Compute min-max shaded ranges across ensembles (dim 1)
min1 = min(er1_all, [], 1); max1 = max(er1_all, [], 1);
min2 = min(er2_all, [], 1); max2 = max(er2_all, [], 1);
min3 = min(er3_all, [], 1); max3 = max(er3_all, [], 1);
% Colors
c_dmd = [0.0000, 0.4470, 0.7410];
c_ked = [0.8500, 0.3250, 0.0980];
c_dual = [0.9290, 0.6940, 0.1250];
% --- 2. Main Plot (All 3 Methods) ---
fill(x_poly, [max2, fliplr(min2)], c_dmd, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(x_poly, [max3, fliplr(min3)], c_ked, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(x_poly, [max1, fliplr(min1)], c_dual, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(mean(er2_all, 1), '--', 'Color', c_dmd, 'LineWidth', 2, 'DisplayName', 'DMD (mean)');
plot(mean(er3_all, 1), '-', 'Color', c_ked, 'LineWidth', 2, 'DisplayName', 'KeDMD (mean)');
plot(mean(er1_all, 1), '-.', 'Color', c_dual, 'LineWidth', 2, 'DisplayName', 'DualKoop-RKHS (mean)');
% Main Axes Formatting
set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 25);
grid on; box on;
title('Relative Forecast Errors Comparison', 'Interpreter', 'latex', 'FontSize', 18);
xlabel('Lead Time (Snapshots)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Relative Forecast Error', 'Interpreter', 'latex', 'FontSize', 18);
legend('Interpreter', 'latex', 'FontSize', 16, 'Location', 'northwest');
% % --- 3. Inset Zoom Plot (KeDMD vs DualKoop-RKHS) ---
% % Position: [left, bottom, width, height]
% axes('Position', [0.50, 0.22, 0.38, 0.38]);
% hold on;
% % Plot shaded areas for KeDMD and DualKoop
% fill(x_poly, [max3, fliplr(min3)], c_ked, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
% fill(x_poly, [max1, fliplr(min1)], c_dual, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
% % Plot mean lines for KeDMD and DualKoop
% plot(mean(er3_all, 1), '-', 'Color', c_ked, 'LineWidth', 2);
% plot(mean(er1_all, 1), '-.', 'Color', c_dual, 'LineWidth', 2);
% % Inset Formatting
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
% grid on; box on;
% title('\textbf{Zoomed: KeDMD vs DualKoop}', 'Interpreter', 'latex', 'FontSize', 11);
% % Automatically adjust inset Y-limits based on KeDMD and DualKoop values
% y_min_zoom = min([min1, min3]);
% y_max_zoom = max([max1, max3]);
% xlim([1, N_steps]);
% ylim([y_min_zoom * 0.9, y_max_zoom * 1.1]);
% --- 4. Export Vector Graphics PDF ---
exportgraphics(gcf, 'ensemble_shallow_weather_error_1.pdf', ...
'ContentType', 'vector', 'BackgroundColor', 'none');

%% Bar plot: Relatev error
% --- 1. Compute overall RMSE across ALL time steps and ALL ensembles ---
% (er(:).^2 squares all elements, mean computes average, sqrt gives RMSE)
rmse_DMD      = sqrt(mean(er2_all(:).^2));
rmse_kEDMD    = sqrt(mean(er3_all(:).^2));
rmse_SpecRKHS = sqrt(mean(er1_all(:).^2));

rmse_values = [rmse_DMD; rmse_kEDMD; rmse_SpecRKHS];

% --- 2. Color definitions (Matching your line plot) ---
c_dmd   = [0.0000 0.4470 0.7410]; % Blue
c_kedmd = [0.8500 0.3250 0.0980]; % Terracotta / Red
c_spec  = [0.9290 0.6940 0.1250]; % Orange

% --- 3. Plot Bar Chart ---
figure('Color', 'w');
b = bar(1:3, rmse_values, 'FaceColor', 'flat', 'BarWidth', 0.5);
b.CData(1,:) = c_dmd;
b.CData(2,:) = c_kedmd;
b.CData(3,:) = c_spec;

% --- 4. Log Scale & Formatting ---
set(gca, 'YScale', 'log');
set(gca, 'XTick', 1:3, 'XTickLabel', {'DMD', 'KeDMD', 'DualKoop-RKHS'});
grid on; box on;
ylabel('Overall RMSE', 'interpreter', 'latex', 'fontsize', 18);
title('Overall Forecast Error Comparison (RMSE)', 'interpreter', 'latex', 'fontsize', 18);
set(gca, 'FontSize', 14);

% Export vector graphics PDF
exportgraphics(gcf, 'ensemble_overall_rmse_bar.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%% ==============Bar ZOOM-------------
figure('Color', 'w', 'Position', [100, 100, 800, 600]);

% --- 1. Main Bar Plot ---
b = bar(1:3, rmse_values, 'FaceColor', 'flat', 'BarWidth', 0.5);
b.CData(1,:) = c_dmd;
b.CData(2,:) = c_kedmd;
b.CData(3,:) = c_spec;

% Apply LaTeX interpreter to X-tick labels and axes
set(gca, 'TickLabelInterpreter', 'latex', ...
    'XTick', 1:3, ...
    'XTickLabel', {'\textrm{DMD}', '\textrm{KeDMD}', '\textrm{DualKoop-RKHS}'});
grid on; box on;
ylabel('Overall RMSE', 'Interpreter', 'latex', 'FontSize', 18);
title('Overall Forecast Error Comparison (RMSE)', 'Interpreter', 'latex', 'FontSize', 18);

% Font size 20 for the main plot axis numbers
set(gca, 'FontSize', 20);

% --- 2. Inset Plot (Narrower width & Boosted zoom) ---
% Position: [left, bottom, width, height] -> Width reduced to 0.22
axes('Position', [0.63, 0.32, 0.22, 0.38]); 
b_inset = bar(3, rmse_values(3), 'FaceColor', 'flat', 'BarWidth', 0.4);
b_inset.CData(1,:) = c_spec;

% Apply LaTeX interpreter to inset plot
set(gca, 'TickLabelInterpreter', 'latex', ...
    'XTick', 3, ...
    'XTickLabel', {'\textrm{DualKoop-RKHS}'});

% Adjust ylim to boost the zoom effect relative to the value magnitude
ylim([rmse_values(3)*0.95, rmse_values(3)*1.05]); 
title('\textbf{Zoomed (DualKoop)}', 'Interpreter', 'latex', 'FontSize', 11);
grid on; box on;

% Use a smaller font size (e.g., 10) for the inset to prevent clipping and text overlap
set(gca, 'FontSize', 13);

% --- 3. Export Vector Graphics PDF ---
exportgraphics(gcf, 'ensemble_overall_rmse_bar_1.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none');

%% =====================================================================
%  USE SpecRKHS-Obs FOR A DOMAIN-MEAN OBSERVABLE
%  (generalizes sealevels.m's "mean sea level" section: instead of the
%  scalar mean height, use the domain-mean of the state field, forecast
%  member-by-member.)
%  =====================================================================
mean_obs_kmd_all   = zeros(steps+1, M);
mean_obs_kedmd_all = zeros(steps+1, M);
mean_obs_exact_all = zeros(steps+1, M);
obs_exact_all      = zeros(steps+1, M);   % exact kernel-space observable
Kx0_er_all         = zeros(steps+1, M);   % relative error in K_{x0} prediction
 
% Domain-mean observable evaluated at the pooled training points X.
% Must have exactly n_train entries (one per column of X, matching G's
% dimension) -- this is fixed across members, so compute it once here,
% mirroring sealevels.m's mean_sea_level_old (indexed over x's columns).
mean_state_train = mean(X, 1).';          % n_train x 1
mode_res = (G * F_res) \ mean_state_train;  % fixed, independent of member
 
for m = 1:M
    x0_m = Xa_3D(:, maxx, m);
 
    % --- recompute per-member SpecRKHS-Obs coefficients at x0_m ---
    Kx0_vals = zeros(n_train, 1);
    for i = 1:n_train
        Kx0_vals(i) = ker(x0_m, X(:, i));
    end
    coefs_res = ((G * F_res) \ Kx0_vals).';
 
    mean_obs_kmd_all(:, m) = real((coefs_res .* (Lambda_res.^(0:steps)).') * ...
        (F_res.' * G * F_res * mode_res)).';
 
    % --- same via kEDMD ---
    G_start = zeros(1, n_train);
    for i = 1:n_train
        G_start(i) = ker(x0_m, X(:, i));
    end
    % Target must match [X, x0_m] used in the main kEDMD section above:
    % n_train pooled columns plus this member's own x0_m appended.
    mean_state_train_full = [mean_state_train; mean(x0_m)];   % (n_train+1) x 1
    mode_full_obs = (([G; G_start] * W) \ mean_state_train_full).';
    psi0_full = G_start * W;
    mean_obs_kedmd_all(:, m) = real(transpose(transpose(psi0_full) .* ...
        (conj(Lambda).^(0:steps))) * mode_full_obs.')';
 
    % --- exact domain-mean over the forecast window (ground truth) ---
    mean_obs_exact_all(:, m) = mean(Ya_3D(:, maxx-1 + (1:steps+1), m), 1).';
 
    % --- exact values of the observable in kernel space (G_future) ---
    G_future = zeros(steps+1, n_train);
    for i = 1:steps+1
        for j = 1:n_train
            G_future(i, j) = ker(Ya_3D(:, maxx-1+i, m), X(:, j));
        end
    end
    obs_exact_all(:, m) = G_future * F_res * mode_res;
 
    % --- error in approximating K_{x0} itself ---
    Kx0_future_vals = zeros(steps+1, 1);
    for i = 0:steps
        Kx0_future_vals(i+1) = ker(x0_m, Ya_3D(:, maxx-1+i+1, m));
    end
    Kx0_predict = real((coefs_res .* (Lambda_res.^(0:steps)).') * F_res.' * Kx0_vals);
    Kx0_er_all(:, m) = abs(Kx0_future_vals - Kx0_predict).^2 ./ abs(Kx0_future_vals).^2;
end
 

%% =====================================================================
%  PLOT COMPARISON OF DOMAIN-MEAN OBSERVABLE PREDICTIONS
%  (averaged over test_members -- a single value if Mt==1)
%  =====================================================================
figure
hold on
p1 = plot(0:steps, log(mean(mean_obs_kedmd_all, 2)), 'linewidth', 3, 'color', [0.8500 0.3250 0.0980]);
p2 = plot(0:steps, log(mean(mean_obs_kmd_all, 2)), '-.','linewidth', 3, 'color', [0.9290 0.6940 0.1250]);
p3 = plot(0:steps, log(mean(mean_obs_exact_all, 2)), '--', 'linewidth', 3, 'color', [0 0.4470 0.7410]);
grid on; box on;
set(gca, 'FontSize', 20);
title('Domain-mean forecast comparison', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Domain-Mean Prediction', 'interpreter', 'latex', 'fontsize', 18)
legend([p3 p1 p2], {'DMD', 'KeDMD', 'DualKoop-RKHS'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'best')
exportgraphics(gcf, 'ensemble_prediction_mean.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')
 
%% relative error: SpecRKHS-Obs observable prediction vs. exact observable
figure
er_obs = mean(abs(obs_exact_all - mean_obs_kmd_all).^2 ./ abs(obs_exact_all).^2, 2);
semilogy(0:steps, (er_obs), 'linewidth', 3, 'color', [0.9290 0.6940 0.1250])
grid on
set(gca, 'FontSize', 20);
title('Domain-mean relative forecast error', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend({'DualKoop-RKHS'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'best')
exportgraphics(gcf, 'ensemble_error_mean.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')
 
%% relative error in predicting the kernel function K_{x0} (mean over M)
figure
semilogy(0:steps, (mean(Kx0_er_all, 2)), 'linewidth', 3, 'color', [0.9290 0.6940 0.1250])
grid on
% Set the font size of the X and Y axis tick labels to 20
set(gca, 'FontSize', 20);
title('Relative forecast error for kernel function', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend({'DualKoop-RKHS'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'best')
exportgraphics(gcf, 'ensemble_prediction_kernel_error.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')
%% =====================================================================
%  Kernel definition (unchanged from your antarctic-based script;
%  sealevels.m uses sigma=1/10000 with a different Matern order -- keep
%  the length-scale that matches your spatial units)
%  =====================================================================
function ker = kernel_matern(x, t)
    sigma = 1/20000;
    r = vecnorm(x - t);
    ker = zeros(1, size(x, 2));
    ker(r>0)  = (sigma*r(r>0)).^(3/2) .* besselk(-3/2, sigma*r(r>0));
    ker(r==0) = sqrt(pi/2);
end


% Main script or main function code goes here


%% Local Function Definition
% function kernel_f = get_kernel(type, X)
%     switch string(type)
%         case "Linear"
%             kernel_f = @(x,y) y'*x;
% 
%         case "Laplacian"
%             d = mean(vecnorm(X - mean(X,2)));
%             if isa(X, 'single')
%                 kernel_f = @(x,y) exp(-pdist2(y', x') / d);
%             else
%                 kernel_f = @(x,y) exp(-sqrt(max(0, -2*real(y'*x) + dot(x,x) + dot(y,y)')) / d);
%             end
% 
%         case "Gaussian"
%             d = mean(vecnorm(X - mean(X,2)));
%             kernel_f = @(x,y) exp(-(-2*real(y'*x) + dot(x,x) + dot(y,y)') / d^2);
% 
%         case "Lorentzian"
%             d = mean(vecnorm(X - mean(X,2)));
%             kernel_f = @(x,y) (1 + (-2*real(y'*x) + dot(x,x) + dot(y,y)') / d^2).^(-1);
% 
%         case "Matern32"
%             d = mean(vecnorm(X - mean(X,2)));
%             kernel_f = @(x,y) (1 + sqrt(3)*sqrt(max(0, -2*real(y'*x) + dot(x,x) + dot(y,y)'))/d) .* ...
%                               exp(-sqrt(3)*sqrt(max(0, -2*real(y'*x) + dot(x,x) + dot(y,y)'))/d);
% 
%         case "Matern52"
%             d = mean(vecnorm(X - mean(X,2)));
%             kernel_f = @(x,y) (1 + sqrt(5)*sqrt(max(0, -2*real(y'*x) + dot(x,x) + dot(y,y)'))/d + ...
%                               (5/3)*(-2*real(y'*x) + dot(x,x) + dot(y,y)')/d^2) .* ...
%                               exp(-sqrt(5)*sqrt(max(0, -2*real(y'*x) + dot(x,x) + dot(y,y)'))/d);
% 
%         otherwise
%             if isnumeric(type)
%                 d = mean(vecnorm(X));
%                 kernel_f = @(x,y) (y'*x/d^2 + 1).^(type);
%             else
%                 error('Unsupported kernel type: %s', string(type));
%             end
%     end
% end