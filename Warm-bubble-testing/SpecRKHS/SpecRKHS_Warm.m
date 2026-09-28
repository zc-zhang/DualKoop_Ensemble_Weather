%% Revisit by means of the SpecRKH.m method
%% https://github.com/GustavConradie1/SpecRKHS/blob/main/examples/antarctic.m

%% Koopman/PF (dual RKHS ResDMD) analysis of your ensemble dataset
% 1:1 structural port of GustavConradie1/SpecRKHS examples/antarctic.m
% Every section of the reference is kept; only the data-loading, kernel,
% and spatial-reshape/mask logic are swapped for your ensemble case.
% clear
% addpath(genpath('./data/your_data'))
% addpath(genpath('./algorithms'))
% addpath(genpath('./colormaps'))
rng(0)


Xa_all = EnVfull_Xa_36;
Ya_all = EnVfull_Ya_36;

p      = size(Xa_all, 1);          % state dimension (reference: d) -> 3880
M      = 10;                       % ensemble members (reference has none: M=1 conceptually)
N_cols = size(Xa_all, 2);
T_snap = N_cols / M;               % time snapshots per ensemble (reference: N)

grid_size = [97, 40];              % [rows(z), cols(y)] -> 97*40 = 3880, matches p
nLAND = (1:p).';                   % identity mask: no land exclusion in your domain

%% =====================================================================
%  BUILD (d,N)-SHAPED SERIES + TRAIN/TEST SPLIT
%  (reference: [d,N]=size(DATA); steps; lag; max; min; x; y)
%  =====================================================================
% Reference works with ONE long trajectory DATA (d x N). Your data is
% ensemble-based, so we first collapse to the ensemble-mean trajectory
% (d x T_snap) to get a single trajectory DATA plays the role of, then
% follow the reference exactly on that trajectory.
Xa_3D = reshape(Xa_all, p, T_snap, M);
Ya_3D = reshape(Ya_all, p, T_snap, M);

DATA_Xa = mean(Xa_3D, 3);          % ensemble-mean input trajectory  (d x T_snap)
DATA_Ya = mean(Ya_3D, 3);          % ensemble-mean output trajectory (d x T_snap)
% For a one-step-ahead system Xa_all(t) -> Ya_all(t), the reference's
% single trajectory DATA(:,t) is best represented here by DATA_Xa; the
% "next state" comes from DATA_Ya at the SAME snapshot index t (not t+1),
% since your (Xa,Ya) pairs are already input/output pairs, unlike the
% reference where y is just x shifted by one column. Keep this in mind:
DATA = DATA_Xa;                    % analogous to reference's DATA
[d, N] = size(DATA);

steps = 6;      % number of held-out future snapshots to forecast (tune)
lag   = 0;      % gap between train/test, analogous to reference's lag=10
max_  = (N - steps - lag);         % renamed from `max` (reserved word in MATLAB)
min_  = 1;

x = DATA_Xa(:, min_:max_-1);       % training input  (reference: x)
y = DATA_Ya(:, min_:max_-1);       % training output (reference: y)
% (reference uses y=DATA(:,min+1:max), i.e. x shifted by one column;
%  here y is the TRUE paired output at the same snapshot, since your
%  Xa/Ya are already (state, next-state) pairs, not a single shifted series)

%% =====================================================================
%  COMPUTE MATRICES   (reference: ker; generate_matrices_kernelized)
%  =====================================================================
% Reference: ker=@(x,t) kernel_matern(x,t);  [G,A,R]=generate_matrices_kernelized(x,y,ker);
% Your kernel_ResDMD call used 'type','Linear' — use the matching kernel here
% so results are directly comparable to your earlier run:
ker = @(x,t) kernel_matern(x,t);
[G, A, R] = generate_matrices_kernelized(x, y, ker);

%% =====================================================================
%  COMPUTE VERIFIED EIGENVALUES
%  (reference: verified_eigenvalues(G,A,R,11))
%  =====================================================================
tol_arg = 11;   % UNVERIFIED: confirm what this argument controls in your
                % local verified_eigenvalues.m (I could not fetch its source)
[Lambda_res, F_res, Lambda, F, res, res_verif, idx, W, W_res] = ...
    verified_eigenvalues(G, A, R, tol_arg);
length(idx)

%% =====================================================================
%  PLOT SPURIOUS AND VERIFIED EIGENVALUES  (verbatim structure)
%  =====================================================================
figure
scatter(angle(Lambda), log(abs(Lambda)), 200, res, '.', 'LineWidth', 1);
hold on
scatter(angle(Lambda_res), log(abs(Lambda_res)), 500, res_verif, '.', 'LineWidth', 1);
box on
clim([0, 0.04])
load('cmap.mat')
colormap(cmap2); colorbar
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
%  COMPUTE PERRON-FROBENIUS (KOOPMAN) MODES AND PLOT
%  (reference: for j=[1 4 6] ... Phi=((G*F)\(x.')).'; imagesc map)
%  =====================================================================
% Reference hand-picks j=[1 4 6] from the ALREADY-VERIFIED set (F_res),
% and reconstructs the spatial mode FRESH for each chosen mode via the
% kernel-space solve — never slicing a pre-built mode matrix. We do the
% same: choose your own indices into the verified set below.
modes_to_show = [1 4 6];   % indices into Lambda_res / F_res (verified set)

for jj = 1:length(modes_to_show)
    j = modes_to_show(jj);
    L_j = Lambda_res(j);
    F_j = F_res(:, j);
    r_j = res_verif(j);

    Phi = ((G * F_j) \ (x.')).';   % reconstruct spatial mode from scratch

    figure
    u = Phi;
    u = real(u * exp(1i * mean(angle(u))));
    u = -u;

    v = zeros(grid_size(1) * grid_size(2), 1) + NaN;
    v(nLAND) = u(:);
    v = reshape(v, grid_size);

    imagesc(v, 'AlphaData', ~isnan(v))
  colormap(ax, brighten(redblueTecplot(21), -0.55));
    colorbar;
    clim([mean(u(:)) - 2*std(u(:)), mean(u(:)) + 2*std(u(:))])
    set(gca, 'Color', [1,1,1]*0.6)
    axis equal
    axis tight
    box on
    set(gca, 'xticklabel', {[]})
    set(gca, 'yticklabel', {[]})
    title(['$\lambda' sprintf('=$ %g+%gi, residual $=$ %f', real(L_j), imag(L_j), r_j)], ...
        'interpreter', 'latex', 'fontsize', 18)
    exportgraphics(gcf, sprintf('ensemble_evals_mode_%d.pdf', j), ...
        'ContentType', 'vector', 'BackgroundColor', 'none', 'resolution', 1000)
end

%% =====================================================================
%  COMPUTE PSEUDOSPECTRA   (verbatim structure)
%  =====================================================================
pts   = 100;
x_pts = linspace(-1.2, 1.2, pts);
y_pts = linspace(-0.02, 1.2, pts/2);
z_pts = kron(x_pts, ones(length(y_pts),1)) + 1i*kron(ones(1,length(x_pts)), y_pts(:));
z_pts = z_pts(:);
res_pspec = pseudospectra(G, A, R, z_pts);

%% =====================================================================
%  PLOT PSEUDOSPECTRAL CONTOURS   (verbatim structure)
%  =====================================================================
res_pspec_rs = reshape(res_pspec, length(y_pts), length(x_pts));
figure
hold on
box on
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
%  SpecRKHS-Obs FORECAST  (reference: x0; Kx0_vals; coefs_res; x_kmd)
%  = the "dual Koopman" forecasting step, using verified eigenpairs only
%  =====================================================================
x0 = DATA_Xa(:, max_);
Kx0_vals = zeros(max_ - min_, 1);
for i = 1:max_ - min_
    Kx0_vals(i) = ker(x0, x(:,i));
end
coefs_res = ((G * F_res) \ Kx0_vals).';
x_kmd = real((conj(coefs_res) .* (conj(Lambda_res).^(1:steps)).') * (F_res' * x(1:d,:).')).';

%% =====================================================================
%  kEDMD BASELINE  (verbatim structure, all eigenpairs, unfiltered)
%  =====================================================================
G_start = zeros(1, max_ - min_);
for i = 1:max_ - min_
    G_start(i) = ker(x0, x(:,i));
end
mode_full = (([G; G_start] * W) \ (DATA_Xa(:, min_:max_).')).';
psi0_full = G_start * W;
x_kedmd = real(transpose(transpose(psi0_full) .* (conj(Lambda).^(1:steps))) * mode_full.')';

%% =====================================================================
%  DMD BASELINE  (verbatim structure)
%  =====================================================================
[U, S, ~] = svd(x, 'econ');
r = rank(S);
U = U(:, 1:r);

PXs = x' * U;
PYs = y' * U;
K = PXs \ PYs;
[W1, LAM, W2] = eig(K, 'vector');
PXr = PXs * W1; PYr = PYs * W1;

c = ([PXr(1,:); PYr]) \ transpose([x y(:,end)]);
x_dmd = real(transpose(transpose(PYr(end,:)) .* (LAM.^(1:steps))) * c)';

%% =====================================================================
%  RELATIVE FORECAST ERROR COMPARISON PLOT  (verbatim structure)
%  =====================================================================
real_data = DATA_Ya(:, N - steps - lag + (1:steps));

er1 = sum(abs(x_kmd   - real_data).^2, 1) ./ sum(abs(real_data).^2, 1);   % SpecRKHS-Obs
er2 = sum(abs(x_dmd   - real_data).^2, 1) ./ sum(abs(real_data).^2, 1);   % DMD
er3 = sum(abs(x_kedmd - real_data).^2, 1) ./ sum(abs(real_data).^2, 1);   % kEDMD

figure
plot(er2, 'linewidth', 2)
hold on
plot(er3, 'linewidth', 2)
plot(er1, 'linewidth', 2)
grid on
title('Relative forecast errors comparison', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend({'DMD','kEDMD','SpecRKHS-Obs'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'best')
exportgraphics(gcf, 'ensemble_forecast_error_delay.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% =====================================================================
%  EXACT MAPS FOR A FEW LEAD-TIME STEPS  (verbatim structure)
%  =====================================================================
for k = 1:min(3, steps)
    figure
    u = real_data(:, k);
    v = zeros(grid_size(1)*grid_size(2), 1) + NaN;
    v(nLAND) = u(:);
    v = reshape(v, grid_size);
    imagesc(v, 'AlphaData', ~isnan(v))
    colormap('parula')
    colorbar
    set(gca, 'Color', [1,1,1]*0.6)
    axis equal
    axis tight
    set(gca, 'xticklabel', {[]})
    set(gca, 'yticklabel', {[]})
    title(sprintf('Exact, step %d', k), 'fontsize', 18, 'interpreter', 'latex')
    exportgraphics(gcf, sprintf('ensemble_exact%d.pdf', k), 'ContentType', 'vector', 'BackgroundColor', 'none')
end

%% =====================================================================
%  PREDICTED MAPS FOR A FEW LEAD-TIME STEPS  (verbatim structure)
%  =====================================================================
for k = 1:min(3, steps)
    figure
    u = real(x_kmd(:, k).');
    v = zeros(grid_size(1)*grid_size(2), 1) + NaN;
    v(nLAND) = u(:);
    v = reshape(v, grid_size');
    imagesc(v, 'AlphaData', ~isnan(v))
    colormap(coolwarm)
    colorbar
    set(gca, 'Color', [1,1,1]*0.6)
    axis equal
    axis tight
    set(gca, 'xticklabel', {[]})
    set(gca, 'yticklabel', {[]})
    title(sprintf('Predicted, step %d', k), 'fontsize', 18, 'interpreter', 'latex')
    exportgraphics(gcf, sprintf('ensemble_predict%d.pdf', k), 'ContentType', 'vector', 'BackgroundColor', 'none')
end

%% =====================================================================
%  BONUS (from your own script, kept): per-ensemble time-series overlay
%  and dominant-mode grid, run AFTER the verified-mode reconstruction
%  above, using KEs/KMs/KEFs from the SAME verified set as x_kmd
%  =====================================================================
% Build primal Koopman objects from the SAME verified eigenpairs used
% above (F_res, Lambda_res), reconstructed fresh, so KMs/KEFs are
% consistent with the residual-verification step -- not the raw 300.
KEs   = Lambda_res;
KEFs  = ((G * F_res) \ (Xa_all.')).';    % verified eigenfunctions evaluated at ALL ensemble points
KMs   = Xa_all * pinv(KEFs.');           % reconstructed fresh, on verified basis only
Ya_pred = real(KMs * diag(KEs) * KEFs.');

indx_spec = 2450;
figure;
plot(Ya_all(indx_spec,:), 'b.-', 'LineWidth', 1.2);
hold on;
plot(real(Ya_pred(indx_spec,:)), 'r-', 'LineWidth', 1.1);

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

%% Dominant-mode spatial grid, sorted, from the VERIFIED set only
highlight_row = 25;
highlight_col = 26;
highlight_index = sub2ind(grid_size, highlight_row, highlight_col);

sorting_criteria = 'Magnitude';
dt = 1;

if strcmp(sorting_criteria, 'RES')
    [~, sort_idx] = sort(res_verif, 'ascend');
elseif strcmp(sorting_criteria, 'Magnitude')
    [~, sort_idx] = sort(abs(KEs), 'descend');
elseif strcmp(sorting_criteria, 'Periodic')
    T_periodic = (2*pi*dt) ./ abs(angle(KEs));
    [~, sort_idx] = sort(T_periodic, 'descend');
end

KEs_sorted  = KEs(sort_idx);
KMs_sorted  = KMs(:, sort_idx);
RES_sorted  = res_verif(sort_idx);

figure;
numRows = 3; numCols = 3;
tiledlayout(numRows, numCols, 'Padding', 'compact', 'TileSpacing', 'compact');
modesToPlot = 1:2:min(17, 2*length(KEs_sorted)-1);

for k = 1:length(modesToPlot)
    mi = modesToPlot(k);
    if mi > length(KEs_sorted); break; end
    KMode_i = abs(KMs_sorted(:, mi));
    mode_i_reshaped = reshape(KMode_i, grid_size');
    ax = nexttile;
    imagesc(mode_i_reshaped);
    colormap(ax, brighten(redblueTecplot(21), -0.55));
    colorbar;
    axis xy; axis equal;
    set(gca, 'FontSize', 10);
    xlabel("y", 'FontSize', 10);
    ylabel("z", 'FontSize', 10);
    lambda_j = KEs_sorted(mi);
    title_str = sprintf('$\\textrm{Mode }%d: \\lambda_{%d} = %.2f %+.2f i$', mi, mi, real(lambda_j), imag(lambda_j));
    title(title_str, 'Interpreter', 'latex', 'FontSize', 10);
    hold on;
    plot(highlight_col, highlight_row, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow');
    hold off;
end
set(gcf, 'Renderer', 'painters');

%% =====================================================================
%  Kernel definitions
%  =====================================================================
% Reference kernel (Matern), kept for reference / comparison if you switch:
function ker = kernel_matern(x,t)
    sigma = 1/20000;
    r = vecnorm(x-t);
    ker = zeros(1, size(x,2));
    ker(r>0)  = (sigma*r(r>0)).^(3/2) .* besselk(-3/2, sigma*r(r>0));
    ker(r==0) = sqrt(pi/2);
end

% % Kernel matching your original kernel_ResDMD 'type','Linear' call:
% function ker = kernel_linear(x,t)
%     ker = x.' * t;
% end

