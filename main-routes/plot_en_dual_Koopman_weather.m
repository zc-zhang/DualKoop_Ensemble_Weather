% plot_en_dual_Koopman_weather.m
%% ============================================================
%% ============================================================
%  ENSEMBLE KOOPMAN ANALYSIS — Full Pipeline
%  Kernel ResDMD on ensemble MSLP data
%  
%  Structure:
%   Sec 1 : Load data
%   Sec 2 : Build reference basis (ONE call to kernel_ResDMD)
%   Sec 3 : Recover UU directly from first-call outputs
%   Sec 4 : Test all ensemble members (no re-call, direct kernel eval)
%   Sec 5 : Plot eigenvalue spectrum
%   Sec 6 : Temporal evolution  phi(x0^m) * lambda^t  per member
%   Sec 7 : Spatial modes in physical space
%   Sec 8 : Residual error bounds + member ranking
%% ============================================================
%   Question (Quantification) Q1 :Does Ensemble mean-field Koopman govern each member directly?
%   Question (Quantification) Q2 : Does ensemble mean-field (dual Koopman) govern the DEVIATION of each member?

%% ============================================================
% 1. Load and prepare data
%% ============================================================
%data_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Ensemble_MSLP';
data_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Ensemble_PREC';
M        = 100;
Xa_ens   = cell(M, 1);
Ya_ens   = cell(M, 1);

%fprintf('Loading %d ensemble members (enMSLP000-099)...\n', M);
fprintf('Loading %d ensemble members (enPREC000-099)...\n', M);
for m = 0:M-1
   % fname       = sprintf('enMSLP%03d.mat', m);
    fname       = sprintf('enPREC%03d.mat', m);
    S_ens       = load(fullfile(data_dir, fname));
    %Xraw        = S_ens.vectorized_MSLP;        % 20860 × 73
     Xraw        = S_ens.vectorized_PREC;        % 20860 × 73
    Xa_ens{m+1} = Xraw(:, 2:end-1);             % 20860 × 71  (x_t)
    Ya_ens{m+1} = Xraw(:, 3:end);               % 20860 × 72  (x_{t+1})
end


%% ================ EnMean MSLP ======================
% fprintf('Loading ensemble mean (enMSLP100)...\n');
% S_mean    = load(fullfile(data_dir, 'enMSLP100.mat'));
% Xraw_mean = S_mean.vectorized_MSLP;             % 20860 × 73

%% ================= En MEan PREC =====================
fprintf('Loading ensemble mean (enPREC100)...\n');
S_mean    = load(fullfile(data_dir, 'enPREC100.mat'));
%Xraw_mean = S_mean.vectorized_MSLP;             % 20860 × 73
Xraw_mean = S_mean.vectorized_PREC;             % 20860 × 73


Xa_mean   = Xraw_mean(:, 2:end-1);              % 20860 × 71
Ya_mean   = Xraw_mean(:, 3:end);                % 20860 × 72

% Geophysical grid
lat_sub        = weatherDat2021AUG_ensemble0.dat_lat;   % 140 × 149
lon_sub        = weatherDat2021AUG_ensemble0.dat_lon;   % 140 × 149
latN           = 140;
lonN           = 149;
shapefile_path = 'D:\Susuki Lab\Testing_Code\data-weather\Data_250401\ne_10m_coastline\ne_10m_coastline.shp';
S_shape        = shaperead(shapefile_path);

%% ============================================================
% 2. Define reference trajectory + parameters
%% ============================================================
% Reference anomaly: member 1 minus ensemble mean
% (change index here if you want a different reference member)
%ref_member = 1;
Xa_ref  = Xa_ens{41,1} - Xa_mean;         % 20860 × 72
Ya_ref  = Ya_ens{41,1} - Ya_mean;         % 20860 × 72
% Xa_ref  = Xa_mean;         % 20860 × 72
% Ya_ref  = Ya_mean;         % 20860 × 72
T_steps = size(Xa_ref, 2);                       % 72

N_dictionary = 300;   % tune: 30 / 50 / 70

%% ============================================================
% 3. Build reference Koopman basis — ONE call only
%% ============================================================
fprintf('\nBuilding reference Koopman basis (N_dict=%d)...\n', N_dictionary);

% Pass Xb = Xa_ref so PSI_x = feature map at training points
[G, K_star, L, PX, PY, PSI_x_ref] = kernel_ResDMD( ...
    Xa_ref, Ya_ref,       ...
    'type', 'Laplacian', ...
    'N',    N_dictionary, ...
    'Xb',   Xa_ref        ...   % PSI_x_ref = G1' * UU = PX
);

% Eigendecomposition of shared dual Koopman operator
[V, D]    = eig(K_star);
Lambda    = diag(D);

% Sort by frequency — slow oscillations first
omega_all = abs(imag(log(Lambda)));
[~, idx]  = sort(omega_all, 'ascend');
Lambda    = Lambda(idx);
V         = V(:, idx);

fprintf('K_star size : %d x %d\n', size(K_star));
fprintf('Top-3 |lambda|: %.4f  %.4f  %.4f\n', ...
    abs(Lambda(1)), abs(Lambda(3)), abs(Lambda(5)));

%% ============================================================
% 4. Test ALL ensemble members — store raw PSI for mode sweep (Q2)
%% ============================================================
fprintf('\nTesting all %d members in fixed reference basis...\n', M);

N_act       = size(K_star, 1);
res_all     = zeros(M, 1);
Phi_ens     = zeros(T_steps, N_act, M);   % real eigenfunction values
PSI_raw_x   = zeros(T_steps, N_act, M);  % raw feature map (before V)
PSI_raw_y   = zeros(T_steps, N_act, M);  % raw feature map Y

for m = 1:M

    %% Q1: =================
     % Xb_m = Xa_ens{m} - Xa_mean;
     % Yb_m = Ya_ens{m} - Ya_mean;

     %% Q2: =================
     Xb_m = Xa_ens{m} ;
     Yb_m = Ya_ens{m};

    [~,~,~,~,~, PSI_xm, PSI_ym] = kernel_ResDMD( ...
        Xa_ref, Ya_ref,       ...
        'type', 'Laplacian',   ...
        'N',    N_dictionary, ...
        'Xb',   Xb_m,         ...
        'Yb',   Yb_m          ...
    );

    res_all(m)        = norm(PSI_ym - PSI_xm * K_star, 'fro') / ...
                        (norm(PSI_xm, 'fro') + 1e-12);
    Phi_ens(:,:,m)    = real(PSI_xm * V);
    PSI_raw_x(:,:,m)  = PSI_xm;    % ← store raw
    PSI_raw_y(:,:,m)  = PSI_ym;    % ← store raw

    if mod(m,10)==0
        fprintf(' Q2 Member %3d/%d   residual = %.4f\n', m, M, res_all(m));
    end
end

%% ============================================================
% 5. Plot: Dual Koopman Eigenvalue Spectrum
%% ============================================================
figure('Position', [100, 100, 700, 600]);
theta = linspace(0, 2*pi, 200);
plot(cos(theta), sin(theta), 'r--', 'LineWidth', 1.5); hold on;
scatter(real(Lambda), imag(Lambda), 80, 'ko', 'filled', 'MarkerEdgeColor', 'k');
hold off;
xlabel('Re(\lambda)', 'FontSize', 12);
ylabel('Im(\lambda)', 'FontSize', 12);
title('Dual Koopman Eigenvalues (shared K^*)', 'FontSize', 14);
%legend('Unit circle', '\lambda_j', 'Location', 'best');
grid on; axis equal;

%% ============================================================
% 6. Temporal evolution: phi(x_0^m) * lambda^t
%    — correct Koopman prediction for each member
%    — compare with direct phi(x_t^m) to show approximation quality
%% ============================================================
eig_plot = 1:2:9;   % which modes to plot
t_vec    = (0:T_steps-1)';

figure('Position', [100, 100, 900, 300*length(eig_plot)]);
for ji = 1:length(eig_plot)
    j = eig_plot(ji);
    subplot(length(eig_plot), 1, ji); hold on;


cmap = jet(M);   % or parula(M), turbo(M) (better than jet sometimes)

for m = 1:M
    phi0_m = Phi_ens(1, j, m);
    phi_pred = real(phi0_m * Lambda(j).^t_vec);

    plot(t_vec, phi_pred, ...
        'Color', cmap(m,:), 'LineWidth', 0.85);
    hold on;
end

colormap(jet);
cb = colorbar;
cb.Label.String = 'Ensemble index m';

    % Reference member prediction (black)
    phi0_ref  = real(PSI_x_ref(1,:) * V(:,j));
    phi_pred_ref = real(phi0_ref * Lambda(j).^t_vec);
    plot(t_vec, phi_pred_ref, 'k-', 'LineWidth', 2.5);

    % Direct evaluation at reference snapshots (dashed red — shows approx error)
    phi_direct_ref = real(PSI_x_ref * V(:,j));
    plot(t_vec, phi_direct_ref, 'r-.', 'LineWidth', 2.5);

    xlabel('Time snapshot t', 'FontSize', 11);
    ylabel(sprintf('\\phi_%d', j), 'FontSize', 11);
    title(sprintf('Mode %d: \\lambda=%.4f+%.4fi, |\\lambda|=%.4f', ...
        j, real(Lambda(j)), imag(Lambda(j)), abs(Lambda(j))), 'FontSize', 11);
    legend('Ensemble \phi(x_0^m)\cdot\lambda^t', ...
           'Reference \phi(x_0)\cdot\lambda^t', ...
           'Direct \phi(x_t) reference', ...
           'Location', 'best', 'FontSize', 9);
    grid on;
end
sgtitle('Dual Koopman Temporal Evolution: \phi(x_0^m)\cdot\lambda^t vs Direct \phi(x_t)', ...
    'FontSize', 13);

%% ============================================================
% 7. Spatial modes in physical space
%    Koopman mode: Xa_ref * PX * V(:,j)  — fixed for all members
%% ============================================================
% eig_spatial = 1:2:17;   % which modes to visualize
% 
% for ji = 1:length(eig_spatial)
%     j = eig_spatial(ji);
% 
%     % Dual Koopman mode in physical space
%     DK_mode   = Xa_ref * (PX * V(:, j));          % 20860 × 1
%     phi_field = reshape(abs(real(DK_mode(1:latN*lonN))), [latN, lonN]);
% 
%     figure('Position', [100, 100, 900, 600]);
% 
%     p = pcolor(lon_sub, lat_sub, phi_field);
%     p.EdgeColor = 'none';
%     shading interp;
%     axis equal tight;
%     set(gca, 'YDir', 'normal');
%     colormap(brighten(redblueTecplot(21), -0.55));
%     cb = colorbar;
%     ylabel(cb, 'Mode Amplitude', 'FontSize', 10);
% 
%     hold on;
%     for k = 1:length(S_shape)
%         plot(S_shape(k).X, S_shape(k).Y, 'k-', 'LineWidth', 1.2);
%     end
%     hold off;
% 
%     xlabel('Longitude (°E)', 'Interpreter', 'tex', 'FontSize', 11);
%     ylabel('Latitude (°N)',  'Interpreter', 'tex', 'FontSize', 11);
%     title(sprintf('Spatial Koopman Mode %d   \\lambda = %.4f + %.4fi   |\\lambda| = %.4f', ...
%         j, real(Lambda(j)), imag(Lambda(j)), abs(Lambda(j))), 'FontSize', 12);
%     xlim([min(lon_sub(:)), max(lon_sub(:))]);
%     ylim([min(lat_sub(:)), max(lat_sub(:))]);
% end


%% Spatial Modes: we use subplot hsere 
eig_spatial = 1:2:17;   % 9 modes (3×3)

figure('Position', [100, 100, 1200, 900]);

tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

for ji = 1:length(eig_spatial)
    j = eig_spatial(ji);

    ax = nexttile;

    % ---- Dual Koopman mode ----
    DK_mode   = Xa_ref * (PX * V(:, j));
    phi_field = reshape(abs(real(DK_mode(1:latN*lonN))), [latN, lonN]);

    % ---- Plot field ----
    p = pcolor(lon_sub, lat_sub, phi_field);
    p.EdgeColor = 'none';
    shading interp;
    axis equal tight;
    set(gca, 'YDir', 'normal');

    colormap(ax, brighten(redblueTecplot(21), -0.55));

    % ---- Individual color scaling ----
    caxis([min(phi_field(:)), max(phi_field(:))]);

    % ---- Proper colorbar (aligned automatically) ----
    cb = colorbar('eastoutside');
    cb.FontSize = 7;   % no label (removed “Amplitude”)

    % ---- Coastline overlay ----
    hold on;
    for k = 1:length(S_shape)
        plot(S_shape(k).X, S_shape(k).Y, 'k-', 'LineWidth', 1.0);
    end
    hold off;

    % ---- Eigenvalue info + period ----
    theta = angle(Lambda(j));
    if abs(theta) > 1e-8
        T = 2*pi / abs(theta);
    else
        T = Inf;
    end

    % ---- One-line title ----
    title(sprintf('V_{%d}: \\lambda=%.3f%+.3fi, |\\lambda|=%.3f, T=%.2f', ...
        j, real(Lambda(j)), imag(Lambda(j)), abs(Lambda(j)), T), ...
        'FontSize', 11);

    % ---- Axis formatting ----
    xlim([min(lon_sub(:)), max(lon_sub(:))]);
    ylim([min(lat_sub(:)), max(lat_sub(:))]);

    if ji > 6
        xlabel('Longitude (°E)', 'FontSize', 8);
    end
    if mod(ji,3)==1
        ylabel('Latitude (°N)', 'FontSize', 8);
    end
end
%% ============================================================
% 8. Residual error bounds + member ranking
%% ============================================================

% --- 8a. Bar chart: residual per member ---
figure('Position', [100, 100, 900, 400]);
bar(1:M, res_all, 'FaceColor', [0.4 0.6 0.9], 'EdgeColor', 'none');
xlabel('Ensemble Member Index', 'FontSize', 12);
ylabel('||PSI_y - PSI_x K^*||_F  /  ||PSI_x||_F', 'FontSize', 11);
title('Residual per Ensemble Member (shared K^*)', 'FontSize', 13);
%yline(mean(res_all), 'r--', 'LineWidth', 2, 'Label', 'Mean');
grid on;

% --- 8b. Ranked bar chart with original member indices ---
[res_sorted, sort_idx] = sort(res_all, 'ascend');
figure('Position', [100, 100, 1200, 500]);
b = bar(res_sorted, 'FaceColor', [0.3 0.75 0.5], 'EdgeColor', 'none');

% Label each bar with original member index
xticks(1:M);
xticklabels(arrayfun(@(i) sprintf('m%d', sort_idx(i)), 1:M, 'UniformOutput', false));
xtickangle(90);

xlabel('Ensemble member (ranked best \rightarrow worst)', 'FontSize', 12);
ylabel('Residual', 'FontSize', 12);
title('Ensemble Members Ranked by Residual', 'FontSize', 13);
grid on;

% Annotate best and worst explicitly
hold on;
text(1,   res_sorted(1),   sprintf(' best\n m%d', sort_idx(1)),   ...
    'FontSize',10, 'Color','b', 'VerticalAlignment','bottom');
text(M,   res_sorted(M),   sprintf(' worst\n m%d', sort_idx(M)),  ...
    'FontSize',10, 'Color','r', 'VerticalAlignment','bottom');
hold off;

fprintf('\nBest  member : m%d  (residual = %.4f)\n', sort_idx(1),   res_sorted(1));
fprintf('Worst member : m%d  (residual = %.4f)\n',   sort_idx(end), res_sorted(end));

% --- 8c. Print summary ---
fprintf('\n========== Residual Summary ==========\n');
fprintf('  Mean : %.4f\n', mean(res_all));
fprintf('  Std  : %.4f\n', std(res_all));
fprintf('  Max  : %.4f  (member %d)\n', max(res_all),  sort_idx(end));
fprintf('  Min  : %.4f  (member %d)\n', min(res_all),  sort_idx(1));
fprintf('\n  Best  5 members: ');  fprintf('%d ', sort_idx(1:5));
fprintf('\n  Worst 5 members: ');  fprintf('%d ', sort_idx(end-4:end));
fprintf('\n=======================================\n');

% --- 8d. Eigenfunction temporal variance across ensemble ---
% Shows which modes have high spread (= uncertain) vs low spread (= robust)
phi_var = squeeze(var(Phi_ens, 0, 1));   % N_act × M  -> var over time dim
% mean variance across members per mode
phi_var_mean = mean(phi_var, 2);         % N_act × 1

figure('Position', [100, 100, 900, 400]);
bar(1:N_act, phi_var_mean, 'FaceColor', [0.85 0.5 0.3], 'EdgeColor', 'none');
xlabel('Eigenfunction (mode) index', 'FontSize', 12);
ylabel('Mean temporal variance across ensemble', 'FontSize', 11);
title('Ensemble Spread in Eigenfunction Coordinates', 'FontSize', 13);
grid on;





%%%%%Additional Box plot

%% ============================================================
% 8c. Boxplot: residual vs number of Koopman modes
%% ============================================================
%% ============================================================
% 8c. Boxplot: residual vs number of Koopman modes
%     As k increases → more modes → better approximation → residual decreases
%     Residual = || PSI_ym_k - PSI_xm_k * K_star_k || 
%     where K_star_k is REFIT (least squares) in k-mode subspace
%% ============================================================
mode_sweep = 1:N_act;
res_modes  = zeros(M, length(mode_sweep));

for ki = 1:length(mode_sweep)
    k   = mode_sweep(ki);
    V_k = V(:, 1:k);                          % top-k eigenvectors

    for m = 1:M
        PSI_xm_k = PSI_raw_x(:,:,m) * V_k;   % T × k
        PSI_ym_k = PSI_raw_y(:,:,m) * V_k;   % T × k

        % REFIT K in k-mode subspace via least squares
        % K_k = argmin || PSI_ym_k - PSI_xm_k * K_k ||
        K_k = PSI_xm_k \ PSI_ym_k;            % k × k  ← refit, not projection

        res_modes(m, ki) = norm(PSI_ym_k - PSI_xm_k * K_k, 'fro') / ...
                           (norm(PSI_xm_k, 'fro') + 1e-12);
    end
end



% % --- Plot boxplot ---
% figure('Position', [100, 100, 1200, 500]);
% boxplot(res_modes, mode_sweep, ...
%     'PlotStyle',   'traditional', ...
%     'MedianStyle', 'line',        ...
%     'Whisker',     1.5,           ...
%     'OutlierSize', 4);
% 
% % Fix x-axis: only show every 5th label to avoid overlap
% ax = gca;
% xticks(5:10:N_act);
% xticklabels(arrayfun(@num2str, 5:5:N_act, 'UniformOutput', false));
% 
% % Reference line: full-mode residual mean
% hold on;
% yline_val = mean(res_all);
% plot([1 N_act], [yline_val yline_val], 'rQuestion 2 — "Does mean-field Koopman govern the DEVIATION of each member?--', 'LineWidth', 2);
% text(N_act*0.7, yline_val*1.05, ...
%     sprintf('Full mode mean = %.3f', yline_val), ...
%     'Color', 'r', 'FontSize', 10);
% hold off;
% 
% xlabel('Number of Koopman modes k', 'FontSize', 12);
% ylabel('Residual over ensemble',     'FontSize', 12);
% title('Residual Error Bound vs Number of Koopman Modes (residual decreases as k increases)', ...
%     'FontSize', 12);
% grid on;
% 
% % --- Median + IQR shaded plot ---
% res_median = median(res_modes, 1);
% res_q1     = quantile(res_modes, 0.25, 1);
% res_q3     = quantile(res_modes, 0.75, 1);
% 
% figure('Position', [100, 100, 1000, 450]);
% fill([mode_sweep, fliplr(mode_sweep)], ...
%      [res_q1,     fliplr(res_q3)],     ...
%      [0.4 0.6 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% hold on;
% plot(mode_sweep, res_median, 'b-', 'LineWidth', 2.5);
% plot([1 N_act], [yline_val yline_val], 'r--', 'LineWidth', 2);
% text(N_act*0.6, yline_val*1.05, ...
%     sprintf('Full basis mean = %.3f', yline_val), ...
%     'Color', 'r', 'FontSize', 10);
% hold off;
% 
% xlabel('Number of Koopman modes k',  'FontSize', 12);
% ylabel('Residual',                   'FontSize', 12);
% title('Median \pm IQR Residual vs Koopman Mode Count', 'FontSize', 13);
% legend('IQR (25th-75th percentile)', 'Median', 'Full basis mean', ...
%     'Location', 'northeast');
% grid on;
% 
% % --- Also plot median + IQR as cleaner summary ---
% res_median = median(res_modes, 1);
% res_q1     = quantile(res_modes, 0.25, 1);
% res_q3     = quantile(res_modes, 0.75, 1);
% 
% figure('Position', [100, 100, 1000, 450]);
% fill([mode_sweep, fliplr(mode_sweep)], ...
%      [res_q1,     fliplr(res_q3)],     ...
%      [0.4 0.6 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
% plot(mode_sweep, res_median, 'b-',  'LineWidth', 2.5);
% %yline(mean(res_all), 'r--', 'LineWidth', 2);
% hold off;
% 
% xlabel('Number of Koopman modes k', 'FontSize', 12);
% ylabel('Residual', 'FontSize', 12);
% title('Median ± IQR Residual vs Koopman Mode Count', 'FontSize', 13);
% legend('IQR (25th–75th)', 'Median', ...
%        sprintf('Full basis mean=%.3f', mean(res_all)), ...
%        'Location', 'best');
% grid on;


%% New part cosnider the 8 parts ---------------

% --- 8a. Bar: residual per member ---
figure('Position', [100 100 900 400]);
bar(1:M, res_all, 'FaceColor', [0.4 0.6 0.9], 'EdgeColor', 'none');
hold on;
plot([0.5 M+0.5],[mean(res_all) mean(res_all)],'r--','LineWidth',2);
text(M*0.75, mean(res_all)*1.03, ...
    sprintf('\\mu=%.3f',mean(res_all)),'Color','r','FontSize',11);
hold off;
xlabel('Ensemble Member Index $i$', ...
    'Interpreter','latex','FontSize',12);
ylabel('$\mathcal{E}^{(i)} = \|\Psi(\delta x^i)K^*_{\bar{x}} - \Psi(\delta y^i)\|_F \,/\, \|\Psi(\delta x^i)\|_F$', ...
    'Interpreter','latex','FontSize',11);
title({'Q2: Does mean-field $K^*$ govern member deviations $\delta x^i = x^i - \bar{x}$?', ...
       '$K^*$ learned from $\bar{x}$,  tested on $\delta x^i$'}, ...
    'Interpreter','latex','FontSize',12);
grid on;

% --- 8b. Ranked bar ---
[res_sorted, sort_idx] = sort(res_all, 'ascend');
figure('Position',[100 100 1200 500]);

% Build color matrix: gray for all, blue for best, red for worst
bar_colors         = repmat([0.7 0.7 0.7], M, 1);
bar_colors(1,   :) = [0.2 0.4 0.8];   % best  → blue
bar_colors(end, :) = [0.8 0.2 0.2];   % worst → red

hold on;
for m = 1:M
    bar(m, res_sorted(m), 'FaceColor', bar_colors(m,:), 'EdgeColor','none');
end
plot([0.5 M+0.5],[mean(res_all)+std(res_all) mean(res_all)+std(res_all)], ...
    'g-.','LineWidth',1.5,'DisplayName','\mu+\sigma threshold');
text(1, res_sorted(1), sprintf(' best\n m%d',sort_idx(1)), ...
    'FontSize',10,'Color',[0.2 0.4 0.8],'VerticalAlignment','bottom');
text(M, res_sorted(M), sprintf(' worst\n m%d',sort_idx(M)), ...
    'FontSize',10,'Color',[0.8 0.2 0.2],'VerticalAlignment','bottom');
hold off;

xticks(1:M);
xticklabels(arrayfun(@(i) sprintf('m%d',sort_idx(i)),1:M,'UniformOutput',false));
xtickangle(90);
xlabel('Ensemble member $i$ ranked $\mathcal{E}^{(i)}$ low $\rightarrow$ high', ...
    'Interpreter','latex','FontSize',12);
ylabel('$\mathcal{E}^{(i)}$','Interpreter','latex','FontSize',13);
title({'Q2: Ensemble members ranked by deviation from mean-field Koopman dynamics', ...
       '{\color{blue}Blue=best},  {\color[rgb]{0.5 0.5 0.5}Gray=others},  {\color{red}Red=worst}'}, ...
    'Interpreter','latex','FontSize',12);
legend('Location','northwest','FontSize',10);
grid on;


% % --- 8c. Boxplot: residual vs number of modes ---
% figure('Position',[100 100 1200 500]);
% boxplot(res_modes, mode_sweep, ...
%     'PlotStyle',   'traditional', ...
%     'MedianStyle', 'line',        ...
%     'Whisker',     1.5,           ...
%     'OutlierSize', 4);
% xticks(5:5:N_act);
% xticklabels(arrayfun(@num2str,5:5:N_act,'UniformOutput',false));
% hold on;
% plot([1 N_act],[1.0 1.0],'k:','LineWidth',1.5);
% text(2,1.02,'No-model baseline','FontSize',9,'Color','k');
% plot([1 N_act],[mean(res_modes(:,1)) mean(res_modes(:,1))],'m--','LineWidth',1.5);
% text(2,mean(res_modes(:,1))*1.03, ...
%     sprintf('k=1 mean=%.3f',mean(res_modes(:,1))),'Color','m','FontSize',9);
% hold off;
% xlabel('Number of Koopman modes $k$','Interpreter','latex','FontSize',12);
% ylabel('$\mathcal{E}^i_k(r) = \|\Psi(\delta x^i)K^*_r - \Psi(\delta y^i)\|_F \,/\, \|\Psi(\delta x^i)\|_F$', ...
%     'Interpreter','latex','FontSize',11);
% title({'Q2: Residual $\mathcal{E}^i_k(r)$ vs number of Koopman modes $r$', ...
%        'More modes $\rightarrow$ richer approximation $\rightarrow$ lower residual'}, ...
%     'Interpreter','latex','FontSize',12);
% grid on;


% --- 8c. Boxplot: residual vs number of modes ---
figure('Position',[100 100 1200 500]);
boxplot(res_modes, mode_sweep, ...
    'PlotStyle',   'traditional', ...
    'MedianStyle', 'line',        ...
    'Whisker',     1.5,           ...
    'OutlierSize', 4);
xticks(5:5:N_act);
xticklabels(arrayfun(@num2str, 5:5:N_act, 'UniformOutput',false));

hold on;
% Meaningful ref 1: full-mode residual mean — this is where curve should END
plot([1 N_act],[mean(res_modes(:,end)) mean(res_modes(:,end))], ...
    'b--','LineWidth',2);
text(N_act*0.5, mean(res_modes(:,end))*0.97, ...
    sprintf('Full mode mean=%.3f', mean(res_modes(:,end))), ...
    'Color','b','FontSize',10);

% Meaningful ref 2: elbow point — optimal k
res_median = median(res_modes, 1);
d_median   = diff(res_median);
[~, elbow] = min(abs(d_median));
plot([elbow elbow], ylim, 'r-', 'LineWidth', 1.5);
text(elbow+0.3, max(res_median)*0.95, ...
    sprintf('Elbow k^*=%d', elbow), ...
    'Color','r','FontSize',10);
hold off;

xlabel('Number of Koopman modes $k$','Interpreter','latex','FontSize',12);
ylabel('$\mathcal{E}^{(i)}(k)$','Interpreter','latex','FontSize',13);
title({'Q2: Residual $\mathcal{E}^{(i)}(k)$ vs Koopman modes $k$', ...
       'Residual decreases as $k\uparrow$,  red line = optimal $k^*$'}, ...
    'Interpreter','latex','FontSize',12);
grid on;

% --- 8d. Median + IQR ---
res_median = median(res_modes,1);
res_q1     = quantile(res_modes,0.25,1);
res_q3     = quantile(res_modes,0.75,1);
d_median   = diff(res_median);
[~, elbow] = min(abs(d_median));

figure('Position',[100 100 1000 450]);
fill([mode_sweep, fliplr(mode_sweep)], ...
     [res_q1,     fliplr(res_q3)],     ...
     [0.4 0.6 0.9],'FaceAlpha',0.3,'EdgeColor','none');
hold on;
plot(mode_sweep, res_median,'b-','LineWidth',2.5);
plot(elbow, res_median(elbow),'r*','MarkerSize',14,'LineWidth',2);
text(elbow+0.5, res_median(elbow), ...
    sprintf('Elbow $k^*=%d$\n$\\mathcal{E}=%.3f$',elbow,res_median(elbow)), ...
    'Interpreter','latex','Color','r','FontSize',10,'VerticalAlignment','bottom');
hold off;
xlabel('Number of Koopman modes $k$','Interpreter','latex','FontSize',12);
ylabel('$\mathcal{E}^{(i)}(k)$','Interpreter','latex','FontSize',13);
title({'Q2: Median $\pm$ IQR of $\mathcal{E}^{(i)}$ — optimal Koopman dimension $k^*$', ...
       '$k^*$ = elbow where residual stops decreasing significantly'}, ...
    'Interpreter','latex','FontSize',12);
legend('IQR (25th--75th percentile)','Median $\mathcal{E}^{(i)}$', ...
    sprintf('Elbow $k^*=%d$',elbow), ...
    'Interpreter','latex','Location','northeast');
grid on;

fprintf('\n===== Q2 Summary =====\n');
fprintf('K* learned from: ensemble mean trajectory\n');
fprintf('Tested on      : member deviations dx^i = x^i - xbar\n');
fprintf('Elbow at k*    : %d modes  (res=%.4f)\n', elbow, res_median(elbow));
fprintf('Best  member   : m%d  E^(i)=%.4f\n', sort_idx(1),   res_sorted(1));
fprintf('Worst member   : m%d  E^(i)=%.4f\n', sort_idx(end), res_sorted(end));
fprintf('Mean E^(i)     : %.4f  Std: %.4f\n', mean(res_all), std(res_all));


%% Accuray plot

%% ============================================================
% Temporal residual: E^(i)(t) per time step
% X = time,  Y = residual,  color = number of modes k
% Shaded = ensemble spread across M members at each time step
%% ============================================================

% mode_select = [10, 20, 50, 65, N_act];
% mode_labels = {'r=10','r=20','r=50','r=65','Full'};
% n_modes     = length(mode_select);
% colors      = lines(n_modes);
% t_vec       = 1:T_steps;
% 
% % Compute per-time-step residual for each mode and each member
% % res_time: T × M × n_modes
% res_time = zeros(T_steps, M, n_modes);
% 
% for ki = 1:n_modes
%     k   = mode_select(ki);
%     V_k = V(:, 1:k);
% 
%     for m = 1:M
%         PSI_xm_k = PSI_raw_x(:,:,m) * V_k;   % T × k
%         PSI_ym_k = PSI_raw_y(:,:,m) * V_k;   % T × k
%         K_k      = PSI_xm_k \ PSI_ym_k;       % k × k refit
% 
%         % Per time step: row-wise residual
%         for t = 1:T_steps
%             num = norm(PSI_ym_k(t,:) - PSI_xm_k(t,:)*K_k, 2);
%             den = norm(PSI_xm_k(t,:), 2) + 1e-12;
%             res_time(t, m, ki) = num / den;
%         end
%     end
% end
% 
% % Mean and std across ensemble members at each time step
%  res_time_mean = squeeze(mean(res_time, 2));   % T × n_modes
% res_time_std  = squeeze(std(res_time,  0,2)); % T × n_modes



mode_select = [10, 20, 50, 65, N_act];
mode_labels = {'r=10','r=20','r=50','r=65','Full'};
n_modes     = length(mode_select);
colors      = lines(n_modes);
t_vec       = 1:T_steps;

res_time = zeros(T_steps, M, n_modes);

for ki = 1:n_modes
    k   = mode_select(ki);
    V_k = V(:, 1:k);

    % Projection onto k-mode subspace — stay in full N-dim space
    P_k = V_k * V_k';             % N×N

    for m = 1:M
        % Full space — do NOT truncate
        PSI_xm = PSI_raw_x(:,:,m);    % T×N
        PSI_ym = PSI_raw_y(:,:,m);    % T×N

        % Per time step: row-wise residual in full space
        for t = 1:T_steps
            pred = PSI_xm(t,:) * P_k * K_star;   % 1×N — k-mode prediction
            targ = PSI_ym(t,:);                   % 1×N — full target
            res_time(t,m,ki) = norm(targ - pred, 2) / ...
                               (norm(targ, 2) + 1e-12);
        end
    end
end

res_time_mean = squeeze(mean(res_time, 2));   % T × n_modes
res_time_std  = squeeze(std(res_time,  0,2)); % T × n_modes
%% --- Plot: temporal residual with shaded ensemble spread ---
figure('Position',[100 100 1100 500]);
hold on;
mode_select = [10, 20, 50, 65, N_act];
mode_labels = {'r=10','r=20','k=50','r=65','Full'};
n_modes     = length(mode_select);
for ki = 1:n_modes
    mu_t  = res_time_mean(:, ki);   % T × 1
    sig_t = res_time_std(:,  ki);   % T × 1

    % Shaded band: mean ± std across ensemble
    fill([t_vec, fliplr(t_vec)],          ...
         [mu_t+sig_t; flipud(mu_t-sig_t)]', ...
         colors(ki,:),                    ...
         'FaceAlpha', 0.15,               ...
         'EdgeColor', 'none',             ...
         'HandleVisibility','off');

    % Mean line
    plot(t_vec, mu_t,              ...
        'Color',     colors(ki,:), ...
        'LineWidth', 2,            ...
        'DisplayName', mode_labels{ki});
end

hold off;
xlabel('Time step $k$',       'Interpreter','latex','FontSize',12);
ylabel('$\mathcal{E}^{(i)}$ mean $\pm$ std over ensemble', ...
    'Interpreter','latex','FontSize',12);
title({'Q2: Temporal residual $\mathcal{E}^{(i)}(t)$ for different Koopman mode counts', ...
       'Shaded = ensemble spread $\pm\sigma$,  Line = ensemble mean'}, ...
    'Interpreter','latex','FontSize',12);
legend('Location','best','FontSize',11);
grid on;

% set(gca, 'YScale', 'log');
% ylabel('$\mathcal{E}^{(i)}(k,t)$ mean $\pm$ std  [log scale]', ...
%     'Interpreter','latex','FontSize',12);


%%-------Log plot--
figure('Position',[100 100 1100 500]);
hold on;

for ki = 1:n_modes
    mu_t  = res_time_mean(:, ki);
    sig_t = res_time_std(:,  ki);

    % Shaded band — clip to positive for log scale
    upper = mu_t + sig_t;
    lower = max(mu_t - sig_t, 1e-6);   % ← prevent log(0) or negative

    fill([t_vec, fliplr(t_vec)],       ...
         [upper; flipud(lower)]',      ...
         colors(ki,:),                 ...
         'FaceAlpha', 0.15,            ...
         'EdgeColor', 'none',          ...
         'HandleVisibility','off');

    plot(t_vec, mu_t,              ...
        'Color',     colors(ki,:), ...
        'LineWidth', 2,            ...
        'DisplayName', mode_labels{ki});
end

% ← log scale: decreasing lines clearly separated
set(gca, 'YScale', 'log');

hold off;
xlabel('Time step $k$','Interpreter','latex','FontSize',12);
ylabel('$\log\,\mathcal{E}^{(i)}$ mean $\pm$ std over ensemble', ...
    'Interpreter','latex','FontSize',12);
title({'Q2: Temporal residual $\mathcal{E}^{(i)}(t)$ — log scale', ...
       'Lower = better,  Shaded = ensemble spread $\pm\sigma$,  Line = mean'}, ...
    'Interpreter','latex','FontSize',12);
legend('Location','best','FontSize',11);
grid on;