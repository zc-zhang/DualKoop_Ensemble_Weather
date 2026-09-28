%% Koopman/PF analysis of ensemble ocean/ice data (adapted from antarctic.m)
%clear
rng(0)

%% =====================================================================
%  LOAD ENSEMBLE DATA (already in memory — no file load needed)
%  =====================================================================
Xa_all = EnVfull_Xa_36;
Ya_all = EnVfull_Ya_36;
DATA_Xa =Xa_all;
DATA_Ya=Ya_all;
p      = size(Xa_all, 1);          % state dimension -> 3880
M      = 10;                       % ensemble members
T_snap = size(Xa_all, 2) / M;      % time snapshots per member
grid_size = [97, 40];              % [rows(z), cols(y)]
nLAND = (1:p).';                   % identity mask: no land exclusion

%% =====================================================================
%  BUILD TRAIN SET (ensemble-aware, reference: x=DATA(:,min:max-1))
%  =====================================================================
steps = 6;      % held-out future snapshots to forecast
lag   = 5;      % gap between train/test window
max_  = T_snap - steps - lag;
min_  = 1;

train_idx = reshape((min_:max_-1)' + (0:M-1)*T_snap, 1, []);

X = Xa_all(:, train_idx);   % p x (n_train*M) -- full ensemble, input
Y = Ya_all(:, train_idx);   % p x (n_train*M) -- full ensemble, output (same-snapshot pair, NOT shifted)

n_train = size(X, 2);        % = (max_-min_)*M -- use this everywhere instead of max_-min_ alone

%% =====================================================================
%  COMPUTE MATRICES   (reference: ker; generate_matrices_kernelized)
%  =====================================================================
ker = @(x,t) kernel_matern(x,t);
[G, A, R] = generate_matrices_kernelized(X, Y, ker);

%% =====================================================================
%  COMPUTE VERIFIED EIGENVALUES  (reference: verified_eigenvalues(G,A,R,11))
%  count mode: always returns exactly `num` lowest-residual eigenpairs
%  =====================================================================
num = 20;   % how many verified eigenpairs to keep -- tune as needed
[Lambda_res, F_res, Lambda, F, res, res_verif, idx, W, W_res] = ...
    verified_eigenvalues(G, A, R, num);   % r omitted -> defaults to full rank
length(idx)

%% =====================================================================
%  PLOT 1: SPURIOUS AND VERIFIED EIGENVALUES  (reference structure)
%  =====================================================================
figure
scatter(angle(Lambda), log(abs(Lambda)), 200, res, '.', 'LineWidth', 1);
hold on
scatter(angle(Lambda_res), log(abs(Lambda_res)), 500, res_verif, '.', 'LineWidth', 1);
box on
clim([0, 0.04])
colormap('parula'); colorbar
xlabel('$\mathrm{arg}(\lambda)$', 'interpreter', 'latex', 'fontsize', 18)
ylabel('$\mathrm{log}(|\lambda|)$', 'interpreter', 'latex', 'fontsize', 18)
title(['Residuals for ensemble data', newline], 'interpreter', 'latex', 'fontsize', 18)
ax = gca; ax.FontSize = 18; axis([-pi pi -20*10^(-3) 10^(-3)])
for k = -5:1:5
    plot(k*pi/6*ones(22,1), -20*10^(-3):10^(-3):10^(-3), '--', 'Color', 'black')
end
xticks([-pi -5*pi/6 -4*pi/6 -3*pi/6 -2*pi/6 -pi/6 0 pi/6 2*pi/6 3*pi/6 4*pi/6 5*pi/6 pi])
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
xticklabels({'$-\pi$','','$-2\pi/3$','','$-\pi/3$','','$0$','','$\pi/3$','','$2\pi/3$','','$\pi$'})
xtickangle(30)
exportgraphics(gcf, 'ensemble_evals_angle.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% =====================================================================
%  PLOT 2: EIGENVALUES ON UNIT CIRCLE (reference-style companion view)
%  =====================================================================
figure;
scatter(real(Lambda), imag(Lambda), 300, res, '.', 'LineWidth', 1);
hold on
plot(cos(0:0.01:2*pi), sin(0:0.01:2*pi), '-k')
axis equal
axis([-1.15, 1.15, -1.15, 1.15])
clim([0, 1])
colormap(jet); colorbar
xlabel('$\mathrm{Re}(\lambda)$', 'interpreter', 'latex', 'fontsize', 18)
ylabel('$\mathrm{Im}(\lambda)$', 'interpreter', 'latex', 'fontsize', 18)
ax = gca; ax.FontSize = 18; box on;
xtickangle(30)
exportgraphics(gcf, 'ensemble_evals_circle.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% =====================================================================
%  PF/KOOPMAN MODE PLOTS  (reference: for j=[1 4 6] — representative
%  verified indices, bounded by however many you actually have)
%  =====================================================================
n_modes_available = length(Lambda_res);
modesToPlot = unique(round(linspace(1, n_modes_available, min(9, n_modes_available))));

figure;
numRows = 3; numCols = 3;
tiledlayout(numRows, numCols, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:length(modesToPlot)
    j = modesToPlot(k);
    L_j = Lambda_res(j);
    F_j = F_res(:, j);
    r_j = res_verif(j);

    Phi = ((G * F_j) \ (X.')).';   % reference formula, applied per verified mode

    u = Phi;
    u = real(u * exp(1i * mean(angle(u))));
    u = -u;

    v = zeros(grid_size(1) * grid_size(2), 1) + NaN;
    v(nLAND) = u(:);
    v = reshape(v, grid_size);

    ax = nexttile;
    imagesc(data.y, data.z, v, 'AlphaData', ~isnan(v))
    colormap(ax, brighten(redblueTecplot(21), -0.55));
    colorbar;
    clim([mean(u(:)) - 2*std(u(:)), mean(u(:)) + 2*std(u(:))])
    set(gca, 'Color', [1,1,1]*0.6)

    axis xy; axis equal
    xlim([0 2e4]); ylim([0 2e4])
    set(gca, 'FontSize', 10);
    xlabel("y", 'FontSize', 10); ylabel("z", 'FontSize', 10);

    title_str = sprintf('$\\textrm{Mode }%d: \\lambda_{%d} = %.2f %+.2f i$', ...
        j, j, real(L_j), imag(L_j));
    title(title_str, 'Interpreter', 'latex', 'FontSize', 10);

    highlight_row = 25; highlight_col = 26;
    hold on;
    plot(data.y(highlight_col), data.z(highlight_row), 'ko', ...
        'MarkerSize', 10, 'MarkerFaceColor', 'yellow');
    hold off;
end
set(gcf, 'Renderer', 'painters');
exportgraphics(gcf, 'ensemble_modes.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% =====================================================================
%  PSEUDOSPECTRA  (reference structure, unchanged)
%  =====================================================================
pts   = 100;
x_pts = linspace(-1.2, 1.2, pts);
y_pts = linspace(-0.02, 1.2, pts/2);
z_pts = kron(x_pts, ones(length(y_pts),1)) + 1i*kron(ones(1,length(x_pts)), y_pts(:));
z_pts = z_pts(:);
res_pspec = pseudospectra(G, A, R, z_pts);

res_pspec_rs = reshape(res_pspec, length(y_pts), length(x_pts));
figure
hold on; box on
v = (10.^(-10:0.3:0));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)), ...
         reshape(imag(z_pts),length(y_pts),length(x_pts)), ...
         log10(real(res_pspec_rs)), log10(v));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)), ...
         -reshape(imag(z_pts),length(y_pts),length(x_pts)), ...
         log10(real(res_pspec_rs)), log10(v));
cbh = colorbar;
cbh.Ticks = log10(10.^(-2:1:0));
cbh.TickLabels = 10.^(-2:1:0);
clim([-2,0]);
reset(gcf)
set(gca, 'YDir', 'normal')
colormap('parula');
axis equal;
plot(sin(0:0.01:2*pi), cos(0:0.01:2*pi), '--', 'color', 'white');
grid on
title('Pseudospectrum of ensemble data', 'interpreter', 'latex', 'fontsize', 18)
xlabel('$\mathrm{Re}(z)$', 'interpreter', 'latex', 'fontsize', 18)
ylabel('$\mathrm{Im}(z)$', 'interpreter', 'latex', 'fontsize', 18)
ax = gca; ax.FontSize = 18; axis equal tight;
axis([x_pts(1), x_pts(end), -y_pts(end), y_pts(end)])
exportgraphics(gcf, 'ensemble_pspec.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% =====================================================================
%  SpecRKHS-Obs + kEDMD FORECAST — PER MEMBER, no ensemble-mean anywhere
%  =====================================================================
x_dkmd_all  = zeros(p, steps, M);
x_kedmd_all = zeros(p, steps, M);

for m = 1:M
    x0_m = Xa_all(:, (m-1)*T_snap + max_);   % member m's own starting snapshot

    % --- SpecRKHS-Obs ---
    Kx0_vals = zeros(n_train, 1);
    for i = 1:n_train
        Kx0_vals(i) = ker(x0_m, X(:,i));
    end
    coefs_res = ((G * F_res) \ Kx0_vals).';
    x_dkmd_all(:,:,m) = real((conj(coefs_res) .* (conj(Lambda_res).^(1:steps)).') * (F_res' * X.')).';

    % --- kEDMD ---
    G_start = zeros(1, n_train);
    for i = 1:n_train
        G_start(i) = ker(x0_m, X(:,i));
    end
    mode_full = (([G; G_start] * W) \ ([X, x0_m].')).';   % [X, x0_m] matches [G; G_start] rows
    psi0_full = G_start * W;
    x_kedmd_all(:,:,m) = real(transpose(transpose(psi0_full) .* (conj(Lambda).^(1:steps))) * mode_full.')';
end

%% =====================================================================
%  DMD BASELINE  (single fit on pooled X,Y -- reference structure)
%  =====================================================================
[U, S, ~] = svd(X, 'econ');
r = rank(S);
U = U(:, 1:r);
PXs = X' * U;
PYs = Y' * U;
K = PXs \ PYs;
[W1, LAM, W2] = eig(K, 'vector');
PXr = PXs * W1; PYr = PYs * W1;
c = ([PXr(1,:); PYr]) \ transpose([X Y(:,end)]);
x_dmd = real(transpose(transpose(PYr(end,:)) .* (LAM.^(1:steps))) * c)';

%% =====================================================================
%  RELATIVE FORECAST ERROR — per member, vs each member's OWN real_data
%  =====================================================================
er1_all = zeros(M, steps);   % SpecRKHS-Obs
er3_all = zeros(M, steps);   % kEDMD

for m = 1:M
    real_data_m = Ya_all(:, (m-1)*T_snap + (T_snap - steps - lag + (1:steps)));
    er1_all(m,:) = sum(abs(x_dkmd_all(:,:,m)  - real_data_m).^2, 1) ./ sum(abs(real_data_m).^2, 1);
    er3_all(m,:) = sum(abs(x_kedmd_all(:,:,m) - real_data_m).^2, 1) ./ sum(abs(real_data_m).^2, 1);
end

% DMD compared against member 1 (pooled fit -> single trajectory; adjust if you
% build a per-member DMD fit instead)
real_data_1 = Ya_all(:, (T_snap - steps - lag + (1:steps)));
er2 = sum(abs(x_dmd - real_data_1).^2, 1) ./ sum(abs(real_data_1).^2, 1);

figure
hold on
plot(log(er2), 'k--', 'linewidth', 2, 'DisplayName', 'DMD')
for m = 1:M
    plot(log(er1_all(m,:)), 'Color', [0.2 0.4 0.8 0.35], 'LineWidth', 1, 'HandleVisibility','off')
    plot(log(er3_all(m,:)), 'Color', [0.9 0.4 0.2 0.35], 'LineWidth', 1, 'HandleVisibility','off')
end
plot(log(mean(er1_all,1)), 'b', 'LineWidth', 2.5, 'DisplayName', 'SpecRKHS-Obs (mean)')
plot(log(mean(er3_all,1)), 'r', 'LineWidth', 2.5, 'DisplayName', 'kEDMD (mean)')
grid on
title('Relative forecast errors comparison (per member + mean)', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('log Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend('Location','best')
exportgraphics(gcf, 'ensemble_forecast_error_delay.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% =====================================================================
%  EXACT vs PREDICTED MAPS FOR A FEW LEAD-TIME STEPS (member 1 shown)
%  =====================================================================
member_to_plot = 1;
real_data_plot = Ya_all(:, (member_to_plot-1)*T_snap + (T_snap - steps - lag + (1:steps)));

for k = 1:min(3, steps)
    figure
    u = real_data_plot(:, k);
    v = zeros(grid_size(1)*grid_size(2), 1) + NaN;
    v(nLAND) = u(:);
    v = reshape(v, grid_size);
    imagesc(data.y, data.z, v, 'AlphaData', ~isnan(v))
    colormap('parula'); colorbar
    set(gca, 'Color', [1,1,1]*0.6)
    axis xy; axis equal; axis tight
    xlim([0 2e4]); ylim([0 2e4])
    set(gca, 'xticklabel', {[]}); set(gca, 'yticklabel', {[]})
    title(sprintf('Exact, member %d, step %d', member_to_plot, k), 'fontsize', 16, 'interpreter', 'latex')
    exportgraphics(gcf, sprintf('ensemble_exact_m%d_%d.pdf', member_to_plot, k), 'ContentType', 'vector', 'BackgroundColor', 'none')
end

for k = 1:min(3, steps)
    figure
    u = real(x_dkmd_all(:, k, member_to_plot));
    v = zeros(grid_size(1)*grid_size(2), 1) + NaN;
    v(nLAND) = u(:);
    v = reshape(v, grid_size);
    imagesc(data.y, data.z, v, 'AlphaData', ~isnan(v))
    colormap('parula'); colorbar
    set(gca, 'Color', [1,1,1]*0.6)
    axis xy; axis equal; axis tight
    xlim([0 2e4]); ylim([0 2e4])
    set(gca, 'xticklabel', {[]}); set(gca, 'yticklabel', {[]})
    title(sprintf('Predicted, member %d, step %d', member_to_plot, k), 'fontsize', 16, 'interpreter', 'latex')
    exportgraphics(gcf, sprintf('ensemble_predict_m%d_%d.pdf', member_to_plot, k), 'ContentType', 'vector', 'BackgroundColor', 'none')
end

%% =====================================================================
%  BONUS: per-member time-series overlay
%  =====================================================================
% K_all: cross-kernel between training points (X) and every point in Xa_all.
% kernel_matern is NOT linear -- use the explicit double loop, not X.'*Xa_all.
N_cols = size(Xa_all, 2);
K_all = zeros(n_train, N_cols);
for i = 1:n_train
    for j = 1:N_cols
        K_all(i,j) = ker(X(:,i), Xa_all(:,j));
    end
end

KEs  = Lambda_res;
KEFs = ((G * F_res) \ K_all).';
KMs  = Xa_all * pinv(KEFs.');
Ya_pred = real(KMs * diag(KEs) * KEFs.');

indx_spec = 2450;
figure;
hold on;
colors = num2cell(lines(M), 2);
for i = 1:M
    start_pt = (i-1)*T_snap + 1;
    end_pt   = i*T_snap;
    xregion(start_pt, end_pt, 'FaceColor', colors{i}, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    text(start_pt + T_snap/2, max(real(Ya_all(indx_spec,:))) * 0.9, ...
        sprintf('Ens# %d', i), 'HorizontalAlignment', 'center', ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', colors{i}*0.7);
end
p1 = plot(Ya_all(indx_spec,:), 'b.-', 'LineWidth', 1.2, 'DisplayName', 'True');
p2 = plot(real(Ya_pred(indx_spec,:)), 'r-', 'LineWidth', 1.1, 'DisplayName', 'Verified Koopman Pred');
grid on;
xlim([1, size(Ya_all, 2)]);
xlabel('Total Snapshot Index (All Ensembles)', 'FontSize', 12);
ylabel(sprintf('State Value at Index %d', indx_spec), 'FontSize', 12);
title(sprintf('Ensemble Trajectory Comparison (State Variable %d)', indx_spec), 'FontSize', 14);
legend([p1, p2], 'Location', 'best');
ax = gca; ax.FontSize = 12;

%% =====================================================================
%  Kernel definition
%  =====================================================================
function ker = kernel_matern(x,t)
    sigma = 1/20000;
    r = vecnorm(x-t);
    ker = zeros(1, size(x,2));
    ker(r>0)  = (sigma*r(r>0)).^(3/2) .* besselk(-3/2, sigma*r(r>0));
    ker(r==0) = sqrt(pi/2);
end