%% Revisit by means of the SpecRKHS.m method
%clear
rng(0)

%% =====================================================================
%  LOAD DATA
%  =====================================================================
Xa_all = EnVfull_Xa_36;
Ya_all = EnVfull_Ya_36;
p      = size(Xa_all, 1);
M      = 10;
T_snap = size(Xa_all, 2) / M;
grid_size = [97, 40];
nLAND = (1:p).';

%% =====================================================================
%  BUILD TRAIN SET
%  =====================================================================
steps = 5;
lag   = 0;
max_  = T_snap - steps - lag;
min_  = 1;

train_idx = reshape((min_:max_-1)' + (0:M-1)*T_snap, 1, []);

% X = Xa_all(:, train_idx);
% Y = Ya_all(:, train_idx);

% use full raw data as tranining data
X = Xa_all;
Y = Ya_all;

n_train = size(X, 2);

%% =====================================================================
%  COMPUTE MATRICES (A here is the DUAL-convention cross matrix, i.e.
%  A(i,j) = kappa(y_j, x_i), as produced by generate_matrices_kernelized)
%  =====================================================================
ker = @(x,t) kernel_matern(x,t);
[G, A, R] = generate_matrices_kernelized(X, Y, ker);
disp('size(G):'); disp(size(G))

%% =====================================================================
%  COMPUTE VERIFIED EIGENVALUES
%  Recall pairing convention out of verified_eigenvalues(G,A,R,...):
%    F, Lambda        -> DUAL eigenpairs      (U* psi = Lambda psi)
%    W, conj(Lambda)  -> PRIMAL eigenpairs    (U  phi = conj(Lambda) phi)
%  =====================================================================
num = 360;
[Lambda_res, F_res, Lambda, F, res, res_verif, idx, W, W_res] = ...
    verified_eigenvalues(G, A, R, num);
disp('length(idx):'); disp(length(idx))
disp('size(W_res):'); disp(size(W_res))   % FIX: check W_res (primal), not F_res

% Primal eigenvalues, correctly paired with W_res:
Lambda_res_primal = conj(Lambda_res);     % FIX: primal eigenvalues need conj()

%% =====================================================================
%  PLOT SPURIOUS AND VERIFIED EIGENVALUES
%  (Residuals were computed for the DUAL operator; if you want to see
%   this in the PRIMAL eigenvalue plane, plot conj(Lambda) instead.
%   Left as dual-plane plot here since res/res_verif were computed
%   for F/F_res, not W/W_res -- see note below.)
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



%%%
figure;
scatter(real(Lambda),imag(Lambda),300,res_verif,'.','LineWidth',1);
hold on
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
axis equal
axis([-1.15,1.15,-1.15,1.15])
clim([0,1])
load('cmap.mat')
colormap(camp); colorbar
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
%title(sprintf('Residuals ($M=%d$)',M),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18; box on;


%% =====================================================================
%  Spatial Visualization: verified Koopman modes
%  FIX: use W_res (primal) instead of F_res (dual); pair with
%       conj(Lambda_res) instead of Lambda_res.
%  =====================================================================

%% PLEASE Revise it tinto ordered of res_indx ? res_verif?
n_modes_available = length(Lambda_res);
%modesToPlot = unique(round(linspace(1, n_modes_available, min(9, n_modes_available))));
modesToPlot =[1,3,5,7,9,11,13,15,17];
figure;
numRows = 3; numCols = 3;
tiledlayout(numRows, numCols, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:length(modesToPlot)
    j = modesToPlot(k);
    L_j = conj(Lambda_res(j));      % FIX: primal eigenvalue
    F_j = W_res(:, j);              % FIX: primal eigenvector

    KMs_Phi = ((G * F_j) \ (X.')).';

    u = KMs_Phi;
    %u = real(u * exp(1i * mean(angle(u))));   % legitimate: fix arbitrary global phase
    % FIX: removed "u = -u;" -- purely cosmetic sign flip, no mathematical basis

    v = zeros(grid_size(1) * grid_size(2), 1) + NaN;
    v(nLAND) = u(:);
    v = reshape(v, grid_size);

    ax = nexttile;
    imagesc(data.y, data.z, abs(v), 'AlphaData', ~isnan(v))
    colormap(ax, brighten(redblueTecplot(21), -0.55));
    colorbar;
    %clim([mean(u(:)) - 2*std(u(:)), mean(u(:)) + 2*std(u(:))])
    set(gca, 'Color', [1,1,1]*0.6)
    axis xy; axis equal
    xlim([0 2e4]); ylim([0 2e4])
    set(gca, 'FontSize', 10);
    xlabel("y", 'FontSize', 10); ylabel("z", 'FontSize', 10);
    title(sprintf('$\\textrm{Mode }%d: \\lambda_{%d} = %.2f %+.2f i$', j, j, real(L_j), imag(L_j)), ...
        'Interpreter', 'latex', 'FontSize', 10);

    hold on;
    plot(data.y(26), data.z(25), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow');
    hold off;
end
set(gcf, 'Renderer', 'painters');
exportgraphics(gcf, 'ensemble_modes.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% =====================================================================
%  PSEUDOSPECTRA
%  FIX: pass A' so KoopPseudoSpec/pseudospectra evaluates the PRIMAL
%  pseudospectrum (its documented convention is <K psi_j, psi_i> = primal).
%  =====================================================================
pts   = 100;
x_pts = linspace(-1.2, 1.2, pts);
y_pts = linspace(-0.02, 1.2, pts/2);
z_pts = kron(x_pts, ones(length(y_pts),1)) + 1i*kron(ones(1,length(x_pts)), y_pts(:));
z_pts = z_pts(:);
res_pspec = pseudospectra(G, A', R, z_pts);   % FIX: A -> A' (primal)

res_pspec_rs = reshape(res_pspec, length(y_pts), length(x_pts));
figure
hold on; box on
v = (10.^(-10:0.3:0));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)), reshape(imag(z_pts),length(y_pts),length(x_pts)), log10(real(res_pspec_rs)), log10(v));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)), -reshape(imag(z_pts),length(y_pts),length(x_pts)), log10(real(res_pspec_rs)), log10(v));
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




%%====================Reconstruction ===============================
% Dual KMD and KMD for Xa_all (or X)


%% =====================================================================
%  SpecRKHS-Obs + kEDMD FORECAST -- ALL MEMBERS, same shared G/F_res/W
%  FIX: SpecRKHS-Obs branch switched from F_res (dual) to W_res (primal).
%       kEDMD branch was already correct (W + conj(Lambda)) -- unchanged.
%  =====================================================================
x_kmd_all   = zeros(p, steps, M);
x_kedmd_all = zeros(p, steps, M);

for m = 1:M
    x0_m = Xa_all(:, (m-1)*T_snap + max_);

    % --- SpecRKHS-Obs (FIXED: F_res -> W_res) ---
    Kx0_vals = zeros(n_train, 1);
    for i = 1:n_train
        Kx0_vals(i) = ker(x0_m, X(:,i));
    end
    coefs_res = ((G * W_res) \ Kx0_vals).';                                  % FIX
    x_kmd_all(:,:,m) = real((conj(coefs_res) .* (conj(Lambda_res).^(1:steps))) * (W_res' * X.')).';  % FIX

    % --- kEDMD (already correct: primal convention throughout) ---
    G_start = zeros(1, n_train);
    for i = 1:n_train
        G_start(i) = ker(x0_m, X(:,i));
    end
    mode_full = (([G; G_start] * W) \ ([X, x0_m].')).';
    psi0_full = G_start * W;
    x_kedmd_all(:,:,m) = real(transpose(transpose(psi0_full) .* (conj(Lambda).^(1:steps))) * mode_full.')';
end

%% =====================================================================
%  DMD BASELINE (unchanged -- standard EDMD, no primal/dual ambiguity
%  here since it's built directly from data, not from a kernel dictionary)
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
%  RELATIVE FORECAST ERROR
%  =====================================================================
er1_all = zeros(M, steps);
er3_all = zeros(M, steps);

for m = 1:M
    real_data_m = Ya_all(:, (m-1)*T_snap + (T_snap - steps - lag + (1:steps)));
    er1_all(m,:) = sum(abs(x_kmd_all(:,:,m)   - real_data_m).^2, 1) ./ sum(abs(real_data_m).^2, 1);
    er3_all(m,:) = sum(abs(x_kedmd_all(:,:,m) - real_data_m).^2, 1) ./ sum(abs(real_data_m).^2, 1);
end

real_data_1 = Ya_all(:, (T_snap - steps - lag + (1:steps)));
er2 = sum(abs(x_dmd - real_data_1).^2, 1) ./ sum(abs(real_data_1).^2, 1);

figure
hold on
plot(er2, 'k--', 'linewidth', 2, 'DisplayName', 'DMD')
for m = 1:M
    plot(er1_all(m,:), 'Color', [0.2 0.4 0.8 0.3], 'LineWidth', 1, 'HandleVisibility', 'off')
    plot(er3_all(m,:), 'Color', [0.9 0.4 0.2 0.3], 'LineWidth', 1, 'HandleVisibility', 'off')
end
plot(mean(er1_all,1), 'b', 'linewidth', 2, 'DisplayName', 'SpecRKHS-Obs (mean)')
plot(mean(er3_all,1), 'r', 'linewidth', 2, 'DisplayName', 'kEDMD (mean)')
grid on
title('Relative forecast errors comparison', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend('location', 'best')
exportgraphics(gcf, 'ensemble_forecast_error_delay.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% =====================================================================
%  EXACT vs PREDICTED MAPS
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
    u = real(x_kmd_all(:, k, member_to_plot));
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
%  BONUS: full-ensemble reconstruction
%  FIX: F_res -> W_res, Lambda_res -> conj(Lambda_res)
%  =====================================================================
%% ---- kernel matrix: training anchors (X) vs. all query points (Xa_all) ----
% N_cols = size(Xa_all, 2);
% K_all = zeros(n_train, N_cols);
% for i = 1:n_train
%     for j = 1:N_cols
%         K_all(i,j) = ker(X(:,i), Xa_all(:,j));
%     end
% end

%% ---- DUAL reconstruction (Lambda_res, F_res are already dual) ----
% Psi_dual     = K_all.' * F_res;            % N_cols x num, dual eigenfunctions at every query point
% Modes_dual   = Xa_all * pinv(Psi_dual.');  % p x num, dual Koopman modes (standard DMD-style fit)
% Ya_pred_dual = real(Modes_dual * diag(Lambda_res) * Psi_dual.');
% 
% %% ---- PRIMAL reconstruction (need conj(Lambda_res), W_res) ----
% Psi_primal   = K_all.' * W_res;
% Modes_primal = Xa_all * pinv(Psi_primal.');
% Ya_pred      = real(Modes_primal * diag(conj(Lambda_res)) * Psi_primal.');
% 
% indx_spec = 2450;
% figure;
% hold on;
% colors = num2cell(lines(M), 2);
% for i = 1:M
%     start_pt = (i-1)*T_snap + 1;
%     end_pt   = i*T_snap;
%     xregion(start_pt, end_pt, 'FaceColor', colors{i}, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
%     text(start_pt + T_snap/2, max(real(Ya_all(indx_spec,:))) * 0.9, ...
%         sprintf('Ens# %d', i), 'HorizontalAlignment', 'center', ...
%         'FontSize', 10, 'FontWeight', 'bold', 'Color', colors{i}*0.7);
% end

%% ---- kernel matrix: training anchors (X) vs. all query points ----
N_cols = size(Xa_all, 2);
K_all = zeros(n_train, N_cols);
for i = 1:n_train
    for j = 1:N_cols
        K_all(i,j) = ker(X(:,i), Xa_all(:,j));
    end
end

%% ---- DUAL (Lambda_res, F_res already dual) ----
KEFs_dual = K_all.' * F_res;              % N_cols x num
KMs_dual  = Xa_all * pinv(KEFs_dual.');   % p x num
KEs_dual  = Lambda_res;                   % num x 1

Ya_pred_dual = real(KMs_dual * diag(KEs_dual) * KEFs_dual.');

%% ---- PRIMAL (need conj(Lambda_res), W_res) ----
KEFs = K_all.' * W_res;                   % N_cols x num
KMs  = Xa_all * pinv(KEFs.');             % p x num
KEs  = conj(Lambda_res);                  % num x 1  <-- FIX: 

Ya_pred = real(KMs * diag(KEs) * KEFs.');

p1 = plot(Ya_all(indx_spec,:), 'b.-', 'LineWidth', 1.2, 'DisplayName', 'True');
p2 = plot(real(Ya_pred_dual(indx_spec,:)), 'r-', 'LineWidth', 1.1, 'DisplayName', 'Verified Koopman Pred');
%p3 = plot(real(Ya_pred(indx_spec,:)), 'g-o', 'LineWidth', 1.1, 'DisplayName', 'Primal Pred');
grid on;
xlim([1, size(Ya_all, 2)]);
xlabel('Total Snapshot Index (All Ensembles)', 'FontSize', 12);
ylabel(sprintf('State Value at Index %d', indx_spec), 'FontSize', 12);
title(sprintf('Ensemble Trajectory Comparison (State Variable %d)', indx_spec), 'FontSize', 14);
legend([p1, p2, p3], 'Location', 'best');
ax = gca; ax.FontSize = 12;

%% =====================================================================
%  Kernel definitions
%  =====================================================================
function ker = kernel_matern(x,t)
    sigma = 1/20000;
    r = vecnorm(x-t);
    ker = zeros(1, size(x,2));
    ker(r>0)  = (sigma*r(r>0)).^(3/2) .* besselk(-3/2, sigma*r(r>0));
    ker(r==0) = sqrt(pi/2);
end