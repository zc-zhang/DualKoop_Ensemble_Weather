%% ============================================================
% DUAL KOOPMAN: FULL DIAGNOSTIC PIPELINE
% Goals:
%   1. Eigenvalue spectrum
%   2. Dual spatial modes
%   3. Dual temporal dynamics
%   4. Dual amplitude distribution (worst/best detection)
%   5. Reconstruction vs truth
%   6. GIF animation
%% ============================================================
%Ya_raw_3D               = reshape(Ya_all, 20860, 70, 100);
%Ya_raw_ens_mean         = mean(Ya_raw_3D, 3);                 % [20860 x 70]
%Ya_mean=Ya_raw_ens_mean;


rng(1);  % reproducible bootstrap

lat_sub = weatherDat2021AUG_ensemble0.dat_lat;  % 140×149
lon_sub = weatherDat2021AUG_ensemble0.dat_lon;  % 140×149
latN = 140;
lonN = 149;
shapefile_path = 'D:\Susuki Lab\Testing_Code\data-weather\Data_250401\ne_10m_coastline\ne_10m_coastline.shp';
S = shaperead(shapefile_path);


%% Specific Loactaion (Kumamoto)
% 2. Calculate index for Kumamoto (32.803N, 130.708E)
targetLat = 32.803; targetLon = 130.708;
dist = sqrt((lat_sub - targetLat).^2 + (lon_sub - targetLon).^2);
[minDist, Kumamoto_Locatidx] = min(dist(:));
[r, c] = ind2sub(size(lat_sub), Kumamoto_Locatidx);
fprintf('Kumamoto is at Index: %d (Row: %d, Col: %d)\n', Kumamoto_Locatidx, r, c);


%% STEP 1: Run kernel_ResDMD (unchanged)
N_dict = 300;
[G, K, L, PX, PY, PSI_x0, PSI_y0, ~, G1, A1, UU, kernel_f] = kernel_ResDMD(...
    Xa_all, Ya_all, ...
    'type',    'Laplacian', ...
    'N',       N_dict, ...
    'cut_off', 1e-6);
    % 'Xb',      Xa_mean, ...  % Out-of-sample input (Ensemble Mean t=0)
    % 'Yb',      Ya_mean);     % Out-of-sample target (Ensemble Mean t=t_1)


%% STEP 2: Eigendecompose K_staryese
% dual Koopman

[V_coeff, Lamda]   = eig(K);
KEs   = diag(Lamda);

%% KEFs
KEFs = PX * V_coeff;  


 %% residual tetsing
   denominators = sum(abs(V_coeff).^2, 1).'; 
  numerators = real(sum(conj(V_coeff) .* (L * V_coeff), 1)).';
  RES = sqrt(max(0, (numerators ./ denominators) - abs(KEs).^2));

%% KMs and reconstruction
 %   % KMs     = Xa_all * pinv(KEFs).';
    KMs = Xa_all * pinv(KEFs.');   %%  KMs = Xa_all / (KEFs.'); correct
 %   % Ya_pred = KMs * Lambda * KEFs.';
  Ya_pred = KMs * diag(KEs) * KEFs.'; 





  %%Calculate and display tracking error
reconstruction_error = norm(Ya_all - real(Ya_pred), 'fro') / norm(Ya_all, 'fro');
fprintf('=======================================\n');
fprintf('Trajectory Reconstruction Error: %.4f%%\n', reconstruction_error * 100);
fprintf('=======================================\n');


  %%  Dual Koopman case 
% K_dual =K_res';
% [V_dual,Lam_dual] = eig(K_dual);
% 
% lam_dual = diag(Lam_dual);
% 
% %% Initial ensemble coordinates
% 
% Phi_x0 = kernel_f(X0_all,Xa_all)'*UU;
% 
% Phi_x0 = Phi_x0.';
% C0 = V_dual \ Phi_x0;
% 
% Phi_modes = (PX*V_coeff).';
% 
% KMs_dual = Xa_all*pinv(Phi_modes);
% 
% M = size(X0_all,2);   % ensemble size = 100
% T = size(Ya_all,2)/M; % time length
% Ya_pred_dual = zeros(size(Ya_all));
% 
% for t = 0:T-1
% 
%     c_t = (lam_dual.^t).*C0;
% 
%     Yt = KMs_dual*c_t;
% 
%     idx=t+1:T:size(Ya_all,2);
% 
%     Ya_pred_dual(:,idx)=Yt;
% 
% end
% 
% %% Global relative error
% 
% err_naive = norm(Ya_all-Ya_pred,'fro') ...
%           / norm(Ya_all,'fro');
% 
% err_dual = norm(Ya_all-Ya_pred_dual,'fro') ...
%          / norm(Ya_all,'fro');
% 
% 
% figure;
% bar([err_naive,err_dual])
% 
% set(gca,'XTickLabel',{'KeDMD','Dual KMD'})
% ylabel('Relative Frobenius error')
% title('Global reconstruction error')
% grid on


%% --- 2. Automated Dual Koopman Sorting Mechanism ---
% %--------Choose your sorting criteria: 'RES', 'Magnitude', or 'Periodic'
sorting_criteria = 'Magnitude'; 
dt = 1; % <--- Set your true Delta t simulation time step here

if strcmp(sorting_criteria, 'RES')
    [~, sort_idx] = sort(RES, 'ascend');
    plot_title_prefix = 'True Mode (Low RES) ';

elseif strcmp(sorting_criteria, 'Magnitude')
    [~, sort_idx] = sort(abs(KEs), 'descend'); 
    plot_title_prefix = 'Dominant Mode ';

elseif strcmp(sorting_criteria, 'Periodic')
    % Direct mathematical evaluation of 2*pi / omega
    T_periodic = (2 * pi * dt) ./ abs(angle(KEs));

    % Sort directly by raw value (oscillations first, non-oscillating/Inf last)
    [~, sort_idx] = sort(T_periodic, 'descend'); 
    plot_title_prefix = 'Periodic Mode ';
end


%RES_sorted= RES(sort_idx);
%Perfectly synchronized re-ordering of all core arrays
KEs  = KEs(sort_idx);        
KMs  = KMs(:, sort_idx);     
KEFs = KEFs(:, sort_idx);
%V_coeff        = V_coeff(:, sort_idx);

% Extra
%N_modes  = min(50, N_dict);
%N_modes=500;
%KEs   = KEs(1:N_modes);
%V_coeff        = V_coeff(:, 1:N_modes);

fprintf('Top 10 dual Koopman eigenvalues:\n');
for j = 1:10
    fprintf('  lambda_%d = %.4f %+.4fi  |lambda|=%.4f\n', ...
        j, real(KEs(j)), imag(KEs(j)), abs(KEs(j)));
end



%% ============================================================
%% PLOT 1: Eigenvalue spectrum
%% ============================================================
figure('Color','w','Position',[100 100 650 580]);
theta = linspace(0, 2*pi, 300);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1.2); hold on;
scatter(real(KEs), imag(KEs), 70, abs(KEs), 'filled');
colormap(jet);
cb = colorbar; ylabel(cb, '|\lambda_j|', 'FontSize', 11);
clim([0 1]);
hold off;
xlabel('Re(\lambda_j)', 'FontSize', 12);
ylabel('Im(\lambda_j)', 'FontSize', 12);
title('Dual Koopman eigenvalue spectrum', 'FontSize', 13);
axis equal; grid on;

%================================
%% Periodic order to extract the modes 
% angles  = abs(angle(Lambda));
% periods = 2*pi ./ max(angles, 1e-10);
% 
% % 2. Sort indices by period (descending: longest period first)
% [~, idx] = sort(periods, 'descend');
% 
% Lambda   = Lambda(idx);
% V        = V(:, idx);
% N_modes  = min(50, N_dict);
% %N_modes=50;
% Lambda   = Lambda(1:N_modes);
% V        = V(:, 1:N_modes);





%% ============================================================
%% 1. DATA & TIME VECTOR SETUP (70 Hours Starting at 14:00)
%% ============================================================
% 70 hourly steps (dropping first 2 zero columns)
startTime = datetime(2021, 8, 10, 14, 0, 0); 
timeVec   = startTime + hours(0:69); % [1 x 70]

% --- A. RAW DATA PREPARATION ---
% Ensemble #1 (20860 x 70)
Ya_raw_ens1             = Ya_all(:, 1:70);
data_raw_ens1_mean      = mean(Ya_raw_ens1, 1);               % [1 x 70]
data_raw_ens1_Kumamoto  = Ya_raw_ens1(Kumamoto_Locatidx, :);  % [1 x 70]

% 100-Member Ensemble Mean (20860 x 70)
Ya_raw_3D               = reshape(Ya_all, 20860, 70, 100);
Ya_raw_ens_mean         = mean(Ya_raw_3D, 3);                 % [20860 x 70]
data_raw_EM_mean        = mean(Ya_raw_ens_mean, 1);           % [1 x 70]
data_raw_EM_Kumamoto    = Ya_raw_ens_mean(Kumamoto_Locatidx, :);% [1 x 70]


% --- B. RECONSTRUCTED DATA PREPARATION ---
% Enforce physical non-negativity constraint on raw Koopman predictions
Ya_pred_clean           = max(real(Ya_pred), 0);

% Ensemble #1 (20860 x 70)
Ya_pred_ens1            = Ya_pred_clean(:, 1:70);
data_pred_ens1_mean     = mean(Ya_pred_ens1, 1);              % [1 x 70]
data_pred_ens1_Kumamoto = Ya_pred_ens1(Kumamoto_Locatidx, :); % [1 x 70]

% 100-Member Ensemble Mean (20860 x 70)
Ya_pred_3D              = reshape(Ya_pred_clean, 20860, 70, 100);
Ya_pred_ens_mean        = mean(Ya_pred_3D, 3);                % [20860 x 70]
data_pred_EM_mean       = mean(Ya_pred_ens_mean, 1);          % [1 x 70]
data_pred_EM_Kumamoto   = Ya_pred_ens_mean(Kumamoto_Locatidx, :);% [1 x 70]

%-----------------------

% % 100-Member dual Ensemble Mean (20860 x 70)
% Ya_pred_dual_clean           = max(real(Ya_pred_dual), 0);
% Ya_pred_dual_3D              = reshape(Ya_pred_dual_clean, 20860, 70, 100);
% Ya_pred_ens_dual_mean        = mean(Ya_pred_dual_3D, 3);                % [20860 x 70]
% data_pred_dual_EM_mean       = mean(Ya_pred_ens_dual_mean, 1);          % [1 x 70]
% data_pred_dual_EM_Kumamoto   = Ya_pred_ens_dual_mean(Kumamoto_Locatidx, :);% [1 x 70]


% --- C. X-AXIS TICK MARK FORMATTING ---
tickTimes  = startTime : hours(6) : (startTime + hours(69));
tickLabels = cell(numel(tickTimes), 1);
for i = 1:numel(tickTimes)
    tickLabels{i} = sprintf('%s %s', datestr(tickTimes(i), 'HH:MM'), datestr(tickTimes(i), 'mm/dd'));
end


%% ============================================================
%% FIGURE 1: LINE PLOT — ENSEMBLE #1 COMPARISON
%% ============================================================
figure('Color', 'w', 'Name', 'Fig 1: Line Plot - Ens #1', 'Position', [100 100 1100 650]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Subplot 1: Kyushu Spatial Mean
ax1 = nexttile;
plot(timeVec, data_raw_ens1_mean, 'b-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Raw (Ens #1)');
hold on;
plot(timeVec, data_pred_ens1_mean, 'r--s', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Reconstructed (Ens #1)');
hold off; grid on;
ylabel('PREC (mm)', 'FontSize', 12);
title('Ensemble #1: Kyushu Spatial Area Average', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 11);
ax1.XTick = tickTimes;
set(ax1, 'XTickLabel', tickLabels, 'XTickLabelRotation', 25, 'FontSize', 11);

% Subplot 2: Kumamoto Local Point
ax2 = nexttile;
plot(timeVec, data_raw_ens1_Kumamoto, 'g-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Raw (Ens #1)');
hold on;
plot(timeVec, data_pred_ens1_Kumamoto, 'm--s', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Reconstructed (Ens #1)');
hold off; grid on;
xlabel('Date / Time (UTC)', 'FontSize', 12);
ylabel('PREC (mm)', 'FontSize', 12);
title('Ensemble #1: Local Rainfall at Kumamoto Point', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 11);
ax2.XTick = tickTimes;
set(ax2, 'XTickLabel', tickLabels, 'XTickLabelRotation', 25, 'FontSize', 11);


%% ============================================================
%% FIGURE 2: GROUPED BAR PLOT — ENSEMBLE #1 COMPARISON
%% ============================================================
figure('Color', 'w', 'Name', 'Fig 2: Bar Plot - Ens #1', 'Position', [150 150 1100 650]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Subplot 1: Kyushu Spatial Mean (Bar)
ax1 = nexttile;
b1 = bar(timeVec, [data_raw_ens1_mean', data_pred_ens1_mean'], 'grouped');
b1(1).FaceColor = [0.2 0.5 0.8]; % Blue for Raw
b1(2).FaceColor = [0.8 0.3 0.3]; % Red for Reconstructed
grid on; ylabel('PREC (mm)', 'FontSize', 12);
title('Ensemble #1: Kyushu Spatial Area Average (Grouped Bars)', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Raw (Ens #1)', 'Reconstructed (Ens #1)'}, 'Location', 'northwest', 'FontSize', 11);
ax1.XTick = tickTimes;
set(ax1, 'XTickLabel', tickLabels, 'XTickLabelRotation', 25, 'FontSize', 11);

% Subplot 2: Kumamoto Local Point (Bar)
ax2 = nexttile;
b2 = bar(timeVec, [data_raw_ens1_Kumamoto', data_pred_ens1_Kumamoto'], 'grouped');
b2(1).FaceColor = [0.1 0.6 0.4]; % Green for Raw
b2(2).FaceColor = [0.7 0.3 0.7]; % Purple for Reconstructed
grid on; xlabel('Date / Time (UTC)', 'FontSize', 12); ylabel('PREC (mm)', 'FontSize', 12);
title('Ensemble #1: Local Rainfall at Kumamoto Point (Grouped Bars)', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Raw (Ens #1)', 'Reconstructed (Ens #1)'}, 'Location', 'northwest', 'FontSize', 11);
ax2.XTick = tickTimes;
set(ax2, 'XTickLabel', tickLabels, 'XTickLabelRotation', 25, 'FontSize', 11);


%% ============================================================
%% FIGURE 3: LINE PLOT — 100-MEMBER ENSEMBLE MEAN COMPARISON
%% ============================================================
figure('Color', 'w', 'Name', 'Fig 3: Line Plot - Ensemble Mean', 'Position', [200 200 1100 650]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Subplot 1: Kyushu Spatial Mean
ax1 = nexttile;
plot(timeVec, data_raw_EM_mean, 'b-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Raw Ensemble Mean');
hold on;
plot(timeVec, data_pred_EM_mean, 'r--s', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Reconstructed Ensemble Mean');
%plot(timeVec, data_pred_dual_EM_mean, 'r--s', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Reconstructed Ensemble Mean');
hold off; grid on;
ylabel('PREC (mm)', 'FontSize', 12);
title('100-Member Ensemble Mean: Kyushu Spatial Area Average', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 11);
ax1.XTick = tickTimes;
set(ax1, 'XTickLabel', tickLabels, 'XTickLabelRotation', 25, 'FontSize', 11);

% Subplot 2: Kumamoto Local Point
ax2 = nexttile;
plot(timeVec, data_raw_EM_Kumamoto, 'g-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Raw Ensemble Mean');
hold on;
%plot(timeVec, data_pred_dual_EM_Kumamoto, 'm--s', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Reconstructed Ensemble Mean');
plot(timeVec, data_pred_EM_Kumamoto, 'm--s', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Reconstructed Ensemble Mean');
hold off; grid on;
xlabel('Date / Time (UTC)', 'FontSize', 12);
ylabel('PREC (mm)', 'FontSize', 12);
title('100-Member Ensemble Mean: Local Rainfall at Kumamoto Point', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 11);
ax2.XTick = tickTimes;
set(ax2, 'XTickLabel', tickLabels, 'XTickLabelRotation', 25, 'FontSize', 11);


%% ============================================================
%% FIGURE 4: GROUPED BAR PLOT — 100-MEMBER ENSEMBLE MEAN COMPARISON
%% ============================================================
figure('Color', 'w', 'Name', 'Fig 4: Bar Plot - Ensemble Mean', 'Position', [250 250 1100 650]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Subplot 1: Kyushu Spatial Mean (Bar)
ax1 = nexttile;
b1 = bar(timeVec, [data_raw_EM_mean', data_pred_EM_mean'], 'grouped');
b1(1).FaceColor = [0.2 0.5 0.8]; % Blue for Raw
b1(2).FaceColor = [0.8 0.3 0.3]; % Red for Reconstructed
grid on; ylabel('PREC (mm)', 'FontSize', 12);
title('100-Member Ensemble Mean: Kyushu Spatial Area Average (Grouped Bars)', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Raw Ensemble Mean', 'Reconstructed Ensemble Mean'}, 'Location', 'northwest', 'FontSize', 11);
ax1.XTick = tickTimes;
set(ax1, 'XTickLabel', tickLabels, 'XTickLabelRotation', 25, 'FontSize', 11);

% Subplot 2: Kumamoto Local Point (Bar)
ax2 = nexttile;
b2 = bar(timeVec, [data_raw_EM_Kumamoto', data_pred_EM_Kumamoto'], 'grouped');
b2(1).FaceColor = [0.1 0.6 0.4]; % Green for Raw
b2(2).FaceColor = [0.7 0.3 0.7]; % Purple for Reconstructed
grid on; xlabel('Date / Time (UTC)', 'FontSize', 12); ylabel('PREC (mm)', 'FontSize', 12);
title('100-Member Ensemble Mean: Local Rainfall at Kumamoto Point (Grouped Bars)', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Raw Ensemble Mean', 'Reconstructed Ensemble Mean'}, 'Location', 'northwest', 'FontSize', 11);
ax2.XTick = tickTimes;
set(ax2, 'XTickLabel', tickLabels, 'XTickLabelRotation', 25, 'FontSize', 11);



%% ================= Raw Ense vs . Ensemble Mean
%% ============================================================
%% 1. SPATIAL ACCUMULATION PREPARATION
%% ============================================================
% Sum precipitation across all 70 time steps (Total Accumulated Rainfall in mm)
Sum_Ya_raw_ens_mean  = sum(Ya_raw_ens_mean, 2);   % [20860 x 1]
Sum_Ya_pred_ens_mean = sum(Ya_pred_ens_mean, 2);  % [20860 x 1]
%Sum_Ya_pred_dual_ens_mean = sum(Ya_pred_ens_dual_mean, 2);

% Absolute Difference / Error Field
Diff_Accumulated = abs(Sum_Ya_raw_ens_mean - Sum_Ya_pred_ens_mean);

% Reshape to 2D Spatial Grids [latN x lonN]
Map_raw  = reshape(Sum_Ya_raw_ens_mean,  [latN, lonN]);
Map_pred = reshape(Sum_Ya_pred_ens_mean, [latN, lonN]);
%Map_dual_pred= reshape(Sum_Ya_pred_dual_ens_mean, [latN, lonN]);
Map_diff = reshape(Diff_Accumulated,    [latN, lonN]);

% Calculate common color limits for Raw vs Predicted
c_min = min([Map_raw(:); Map_pred(:)]);
c_max = max([Map_raw(:); Map_pred(:)]);


%% ============================================================
%% 2. FIGURE CREATION: 3-PANEL SPATIAL ACCUMULATION MAPS
%% ============================================================
figure('Color', 'w', 'Name', 'Spatial Accumulated Rainfall Comparison', 'Position', [100 100 1400 450]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% -------------------------------------------------------------
% SUBPLOT 1: RAW ENSEMBLE MEAN
% -------------------------------------------------------------
nexttile;
p1 = pcolor(lon_sub, lat_sub, Map_raw);
p1.EdgeColor = 'none'; shading interp; axis equal tight;
set(gca, 'YDir', 'normal', 'FontSize', 10);
clim([c_min, c_max]);
colormap(brighten(redblueTecplot(21),-0.55));
%colormap(gca, parula); % Standard precipitation colormap
cb1 = colorbar; cb1.Label.String = 'Total PREC (mm)';
hold on;
% Draw Coastlines
for kk = 1:length(S)
    plot(S(kk).X, S(kk).Y, 'k-', 'LineWidth', 1.2);
end
% Mark Kumamoto Point
plot(lon_sub(r,c), lat_sub(r,c), 'p', ...
    'MarkerSize', 12, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
hold off;
title('Raw Ensemble Mean (\Sigma PREC)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Lon (\circE)'); ylabel('Lat (\circN)');
xlim([min(lon_sub(:)), max(lon_sub(:))]);
ylim([min(lat_sub(:)), max(lat_sub(:))]);

% -------------------------------------------------------------
% SUBPLOT 2: RECONSTRUCTED ENSEMBLE MEAN
% -------------------------------------------------------------
nexttile;
p2 = pcolor(lon_sub, lat_sub, Map_pred);
%p2 = pcolor(lon_sub, lat_sub, Map_dual_pred);
p2.EdgeColor = 'none'; shading interp; axis equal tight;
set(gca, 'YDir', 'normal', 'FontSize', 10);
clim([c_min, c_max]); % Enforce identical color scale as Raw
%colormap(gca, parula);
colormap(brighten(redblueTecplot(21),-0.55));
cb2 = colorbar; cb2.Label.String = 'Total PREC (mm)';
hold on;
for kk = 1:length(S)
    plot(S(kk).X, S(kk).Y, 'k-', 'LineWidth', 1.2);
end
plot(lon_sub(r,c), lat_sub(r,c), 'p', ...
    'MarkerSize', 12, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
hold off;
title('Reconstructed Ensemble Mean (\Sigma PREC)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Lon (\circE)'); ylabel('Lat (\circN)');
xlim([min(lon_sub(:)), max(lon_sub(:))]);
ylim([min(lat_sub(:)), max(lat_sub(:))]);

% -------------------------------------------------------------
% SUBPLOT 3: ABSOLUTE RECONSTRUCTION ERROR
% -------------------------------------------------------------
nexttile;
p3 = pcolor(lon_sub, lat_sub, Map_diff);
p3.EdgeColor = 'none'; shading interp; axis equal tight;
set(gca, 'YDir', 'normal', 'FontSize', 10);
colormap(brighten(redblueTecplot(21),-0.55));
cb3 = colorbar; cb3.Label.String = 'Absolute Error (mm)';
hold on;
for kk = 1:length(S)
    plot(S(kk).X, S(kk).Y, 'k-', 'LineWidth', 1.2);
end
plot(lon_sub(r,c), lat_sub(r,c), 'p', ...
    'MarkerSize', 12, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
hold off;
title('|Raw - Reconstructed| Error', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Lon (\circE)'); ylabel('Lat (\circN)');
xlim([min(lon_sub(:)), max(lon_sub(:))]);
ylim([min(lat_sub(:)), max(lat_sub(:))]);








%% ============================================================
%% PLOT: Top 9 dual spatial modes order by the amgnitudes
%% ============================================================
figure('Color','w','Position',[100 100 1200 900]);
colormap(brighten(redblueTecplot(21), -0.55));   % set once

p = 1; % Initialize the subplot position counter
for k = 1:2:17
    mode_map = reshape(abs(KMs(:,k)), [latN, lonN]);
    clim_m   = max(abs(mode_map(:)));

    c_max = max(mode_map(:));
    c_min = min(mode_map(:)); % Keep min if there are negative anomalies in your modes

    ang      = abs(angle(KEs(k)));
    period   = 2*pi / max(ang, 1e-10);

    % --- THIS IS THE KEY CHANGE ---
    subplot(3,3,p); 
    % ------------------------------
    p_plot = pcolor(lon_sub, lat_sub, mode_map);
    p_plot.EdgeColor = 'none'; shading interp; axis equal tight;
    set(gca,'YDir','normal');
   clim([c_min, c_max]);
    colorbar;
    hold on;
    for kk = 1:length(S)
        plot(S(kk).X, S(kk).Y, 'k-', 'LineWidth', 1.5);
    end

    %% mark Kumamoto
    plot(lon_sub(r,c), lat_sub(r,c), 'p', ...
        'MarkerSize',12,'MarkerFaceColor','g','MarkerEdgeColor','k');
    hold off;
    title(sprintf('$V_{%d}$: $|\\lambda|=%.3f$, $T=%.1f$', ...
        k, abs(KEs(k)), period), ...
        'Interpreter','latex','FontSize',10);
    xlabel('Lon ^\circE'); ylabel('Lat ^\circN');
    xlim([min(lon_sub(:)), max(lon_sub(:))]);
    ylim([min(lat_sub(:)), max(lat_sub(:))]);

    p = p + 1; % Increment p so the next mode goes into the next subplot slot
end
sgtitle('Dual Koopman spatial modes $V_j$', ...
    'Interpreter','latex','FontSize',14,'FontWeight','bold');


%%  ===============In-sample vs out-of-sample 
%% ========================================================================
%% Nr sweep: two error bounds — in-sample (train) and out-of-sample (mean)
%% ========================================================================
KEs_sorted     = KEs(sort_idx);
V_coeff_sorted = V_coeff(:, sort_idx);

%% ---- STEP 1: Fit KMs and KEFs ONCE, full rank, on training data ----
KEFs_sorted = PX * V_coeff_sorted;                      % training eigenfunctions, full N_dict
pinv_tol    = 1e-2;                                       % check this filters properly for your data
KMs_sorted  = Xa_all * pinv(KEFs_sorted.', pinv_tol);      % fit once, full rank

% Out-of-sample eigenfunctions (fixed, evaluated once on ensemble mean point)
KEFs_mean_sorted = PSI_x0 * V_coeff_sorted;

%% ---- Sanity check: confirm sample-vs-feature orientation ----
fprintf('size(Xa_all)  = [%d, %d]\n', size(Xa_all));
fprintf('size(Ya_all)  = [%d, %d]\n', size(Ya_all));
fprintf('size(Ya_mean) = [%d, %d]\n', size(Ya_mean));
% Ya_all and Xa_all should share the same orientation (features x samples)

%% ========================================================================
%% STEP 2: Nr sweep -- one-shot reconstruction, no refit per r
%% ========================================================================
r_values =[50, 100, 150, 200, 250, 300];
num_r    = length(r_values);

% Determine sample axis automatically: samples = matching dim to Xa_all
n_features_Ya = size(Ya_all, 1);
n_samples_Ya  = size(Ya_all, 2);   % assume features x samples (standard convention)

err_in_sample_mean = zeros(num_r, 1);
err_in_sample_std   = zeros(num_r, 1);
err_out_mean        = zeros(num_r, 1);
err_out_std          = zeros(num_r, 1);

for k = 1:num_r
    r = r_values(k);

    KMs_r  = KMs_sorted(:, 1:r);
    KEs_r  = KEs_sorted(1:r);
    KEFs_r = KEFs_sorted(:, 1:r);
    KEFs_mean_r = KEFs_mean_sorted(:, 1:r);

    % --- In-sample reconstruction: Xa_all -> Ya_all (all training snapshots) ---
    Ya_pred_r = real(KMs_r * diag(KEs_r) * KEFs_r.');   % features x samples

    % Per-snapshot (per-column) relative error -> distribution across samples
    snap_err = vecnorm(Ya_all - Ya_pred_r, 2, 1) ./ vecnorm(Ya_all, 2, 1);  % 1 x n_samples
    err_in_sample_mean(k) = mean(snap_err);
    err_in_sample_std(k)  = std(snap_err);

    % --- Out-of-sample reconstruction: Xa_mean -> Ya_mean ---
    Ya_pred_mean_r = real(KMs_r * diag(KEs_r) * KEFs_mean_r.');  % features x 1 (or few)

    % Per-element (per-feature/spatial-point) relative error -> distribution across features
    elem_err = abs(Ya_mean(:) - Ya_pred_mean_r(:)) ./ (abs(Ya_mean(:)) + eps);
    err_out_mean(k) = mean(elem_err);
    err_out_std(k)  = std(elem_err);
end

disp(table(r_values(:), err_in_sample_mean, err_in_sample_std, err_out_mean, err_out_std, ...
    'VariableNames', {'r','err_in_mean','err_in_std','err_out_mean','err_out_std'}))

%% ========================================================================
%% STEP 3: Two-panel error-bar figure
%% ========================================================================
fig = figure('Units','inches','Position',[1,1,10,4.2],'PaperPositionMode','auto');

subplot(1,2,1);
errorbar(r_values, err_in_sample_mean, err_in_sample_std, '-s', ...
    'Color',[0.20,0.40,0.75], 'LineWidth',1.8, 'MarkerSize',7, ...
    'MarkerFaceColor',[0.20,0.40,0.75], 'CapSize',6);
set(gca,'YScale','log','FontName','Times New Roman','FontSize',10); grid on;
xlabel('Retained Modes (r)','FontWeight','bold');
ylabel('Relative Error (mean \pm std)','FontWeight','bold');
title('In-Sample: Ensemble Snapshot Spread','FontWeight','bold');

subplot(1,2,2);
errorbar(r_values, err_out_mean, err_out_std, '-o', ...
    'Color',[0.85,0.30,0.25], 'LineWidth',1.8, 'MarkerSize',7, ...
    'MarkerFaceColor',[0.85,0.30,0.25], 'CapSize',6);
set(gca,'YScale','log','FontName','Times New Roman','FontSize',10); grid on;
xlabel('Retained Modes (r)','FontWeight','bold');
ylabel('Relative Error (mean \pm std)','FontWeight','bold');
title('Out-of-Sample: Spatial/Feature Spread','FontWeight','bold');

sgtitle('In-Sample vs. Out-of-Sample Error Bounds Across Mode Truncation', ...
    'FontName','Times New Roman','FontSize',12,'FontWeight','bold');

%% ==================Anomonly correctiion ciefficinet 
%% Out-of-Sample Dual Amplitudes (RKHS Perturbation Analysis)
% Phi0: Eigenfunction evaluations at out-of-sample initial states [M x N_modes]
Phi0      = PSI_x0 * V_coeff;               % phi_j(x_0^i)
Phi0_mean = mean(Phi0, 1);                  % phi_j(x_0_mean)

% Alpha_all: Initial dual amplitude perturbations relative to the ensemble mean
Alpha_all = Phi0 - Phi0_mean;               % [M x N_modes]: delta_alpha_j^i

% Get dynamic dimensions
[M_members, N_modes] = size(Phi0);          % M = number of ensemble members in X0_all

%% STEP 4: Physical Spatial Modes
% Project dual eigenvector coefficients back into physical state space
Modes_space=KMs;
%Modes_space = Xa_all * (PX * V_coeff);      % [p x N_modes]

% Optional normalization for spatial mode visualization
Modes_norm  = Modes_space ./ max(abs(Modes_space), [], 1);

%% STEP 5: Temporal Evolution per Mode Across All Ensemble Members
K_forecast   = 35;                          % Forecast time horizon
t_vec        = (0:K_forecast)';             % [K_forecast+1 x 1]

% Tensor shape: [Time_Steps x Ensemble_Members x Modes]
temporal_all = zeros(K_forecast + 1, M_members, N_modes);

for j = 1:N_modes
    lam_k = (KEs(j)) .^ t_vec;              % Eigenvalue time progression [K_forecast+1 x 1]
    
    % Temporal trajectory for mode j across all M members
    temporal_all(:, :, j) = real(lam_k * Phi0(:, j)');  % [(K_forecast+1) x M_members]
end

%% STEP 6: RKHS Anomaly Score (Hilbert Space Risk Metric)
% ||delta_kappa_k^i||_H^2 = sum_j |lambda_j|^{2k} * |alpha_j^i|^2
anomaly_score = zeros(M_members, K_forecast + 1);  % [M_members x (K_forecast+1)]

for k = 0:K_forecast
    weights = (abs(KEs)) .^ (2 * k);        % Mode energy weights at time step k [N_modes x 1]
    
    % Squared RKHS perturbation norm for each member i at step k
    anomaly_score(:, k + 1) = (abs(Alpha_all) .^ 2) * weights;
end

%% STEP 7: Reconstruct Physical Forecasts (Optional - For Plotting Maps)
% Reconstruct full physical state Y_pred for all members over time: [p x Time_Steps x Members]
p_points = size(Modes_space, 1);
Y_pred_ensemble = zeros(p_points, K_forecast + 1, M_members);

for t_idx = 1:(K_forecast + 1)
    % Mode sum at time step t_idx for all members
    % Modes_space: [p x N_modes], temporal: [M_members x N_modes]
    temp_t = squeeze(temporal_all(t_idx, :, :)); % [M_members x N_modes]
    Y_pred_ensemble(:, t_idx, :) = Modes_space * temp_t'; 
end

%% STEP 8: Detect & Plot Worst-Case / Highest-Risk Member
figure('Color', 'w', 'Position', [100, 100, 800, 450]);

% Plot all ensemble perturbation curves in gray
plot(0:K_forecast, anomaly_score', 'Color', [0.7, 0.7, 0.7, 0.4], 'LineWidth', 1);
hold on;

% Highlight Ensemble Mean Anomaly (0 reference baseline)
[~, worst_idx] = max(anomaly_score(:, end)); % Find member with highest error at t_end
plot(0:K_forecast, anomaly_score(worst_idx, :), 'r-', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Worst-Case Member (#%d)', worst_idx));

grid on; box on;
xlabel('Forecast Time Step (k)', 'FontSize', 11);
ylabel('RKHS Anomaly Score  ||\Delta\kappa_k^i||_{\mathcal{H}}^2', 'FontSize', 11);
title('Dual Koopman Out-of-Sample Ensemble Perturbation Growth', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10);


%% ============================================================
%% PLOT 3 REVISED: Temporal dynamics — full superposition at Kumamoto
%% compared to truth Xa_ens
%% ============================================================
M=100;
% Grid point index for Kumamoto (already defined as r,c in your code)
pt_idx = sub2ind([latN, lonN], r, c);   % linear index into p-vector

% --- Dual KMD reconstruction time series at Kumamoto ---
% recon_ts(k, i) = sum_j Re(lambda_bar_j^k * delta_alpha_j^i * V_j(pt_idx))
recon_ts = zeros(K_forecast+1, M);
for k_t = 0:K_forecast
    %lam_k = conj(Lambda).^k_t;          % N_modes × 1
     lam_k = (KEs).^k_t;        
    % full superposition: p × M, then extract Kumamoto row
    recon_ts(k_t+1, :) = real(...
        Modes_space(pt_idx,:) * diag(lam_k) * Alpha_all');
    % shape: 1 × M
end

% --- Truth time series at Kumamoto ---
truth_ts = zeros(K_forecast+1, M);
for m = 1:M
    for k_t = 0:K_forecast
        truth_ts(k_t+1, m) = ...
            Xa_ens{m}(pt_idx, k_t+1) - Xa_mean(pt_idx, k_t+1);
    end
end

% --- Identify worst/best by RKHS anomaly score at final time ---
[~, worst_m] = max(anomaly_score(:, end));
[~, best_m]  = min(anomaly_score(:, end));

% --- Plot ---
figure('Color','w','Position',[100 100 1200 500]);

% Panel 1: Dual KMD reconstruction
subplot(1,2,1);
hold on;
for m = 1:M
    plot(t_vec, recon_ts(:,m), ...
        'Color',[0.6 0.6 0.6 0.25],'LineWidth',0.5);
end
plot(t_vec, recon_ts(:,worst_m), 'r-','LineWidth',2.5,...
    'DisplayName',sprintf('Worst $m=%d$',worst_m));
plot(t_vec, recon_ts(:,best_m),  'b-','LineWidth',2.5,...
    'DisplayName',sprintf('Best $m=%d$',best_m));
plot(t_vec, mean(recon_ts,2), 'k-','LineWidth',2,...
    'DisplayName','Ensemble mean');
t_std = std(recon_ts,[],2);
t_mean = mean(recon_ts,2);
fill([t_vec; flipud(t_vec)], ...
    [t_mean+t_std; flipud(t_mean-t_std)], ...
    'k','FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
hold off;
xlabel('Forecast step $k$','Interpreter','latex','FontSize',12);
ylabel('PREC deviation (Pa)','FontSize',12);
title('Dual KMD: $\sum_j \mathrm{Re}(\bar\lambda_j^k\,\delta\alpha_j^i\,V_j)$ at Kumamoto',...
    'Interpreter','latex','FontSize',11);
legend('Interpreter','latex','Location','best','FontSize',9);
grid on;

% Panel 2: Truth from Xa_ens
subplot(1,2,2);
hold on;
for m = 1:M
    plot(t_vec, truth_ts(:,m), ...
        'Color',[0.6 0.6 0.6 0.25],'LineWidth',0.5);
end
plot(t_vec, truth_ts(:,worst_m), 'r-','LineWidth',2.5,...
    'DisplayName',sprintf('Worst $m=%d$',worst_m));
plot(t_vec, truth_ts(:,best_m),  'b-','LineWidth',2.5,...
    'DisplayName',sprintf('Best $m=%d$',best_m));
plot(t_vec, mean(truth_ts,2), 'k-','LineWidth',2,...
    'DisplayName','Ensemble mean');
t_std2  = std(truth_ts,[],2);
t_mean2 = mean(truth_ts,2);
fill([t_vec; flipud(t_vec)], ...
    [t_mean2+t_std2; flipud(t_mean2-t_std2)], ...
    'k','FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
hold off;
xlabel('Forecast step $k$','Interpreter','latex','FontSize',12);
ylabel('RREC deviation (mm) ','FontSize',12);
title('Truth: $y_k^i - \bar{y}_k$ at Kumamoto',...
    'Interpreter','latex','FontSize',11);
legend('Interpreter','latex','Location','best','FontSize',9);
grid on;

sgtitle(sprintf(['Dual KMD superposition vs truth at Kumamoto'...
    ' — worst/best by RKHS score, $N_{\\mathrm{modes}}=%d$'],...
    N_modes),...
    'Interpreter','latex','FontSize',13,'FontWeight','bold');
%% ============================================================
%% PLOT 4: Dual amplitude distribution
%% ============================================================
figure('Color','w','Position',[100 100 1400 500]);
% Panel 1: |alpha_j^i| heatmap
subplot(1,3,1);
imagesc(1:N_modes, 1:M, abs(Alpha_all));
colorbar; colormap(brighten(redblueTecplot(21),-0.55));
xlabel('Mode $j$','Interpreter','latex','FontSize',12);
ylabel('Ensemble member $i$','Interpreter','latex','FontSize',12);
title('$|\delta\alpha_j^i|$','Interpreter','latex','FontSize',13);
% mark worst/best
[~, worst_m] = max(anomaly_score(:,end));
[~, best_m]  = min(anomaly_score(:,end));
hold on;
yline(worst_m,'r-','LineWidth',1.5,'DisplayName','Worst');
yline(best_m, 'b-','LineWidth',1.5,'DisplayName','Best');
hold off;

% Panel 2: Re(alpha_j^i) signed
subplot(1,3,2);
clim_a = max(abs(real(Alpha_all(:))));
imagesc(1:N_modes, 1:M, real(Alpha_all));
colorbar;
clim([-clim_a, clim_a]);
xlabel('Mode $j$','Interpreter','latex','FontSize',12);
ylabel('Ensemble member $i$','Interpreter','latex','FontSize',12);
title('$\mathrm{Re}(\delta\alpha_j^i)$','Interpreter','latex','FontSize',13);
hold on;
yline(worst_m,'r-','LineWidth',1.5);
yline(best_m, 'b-','LineWidth',1.5);
hold off;

% Panel 3: RKHS anomaly score over time
subplot(1,3,3);
hold on;
for m = 1:M
    plot(t_vec, anomaly_score(m,:), ...
        'Color',[0.7 0.7 0.7 0.2],'LineWidth',0.5);
end
plot(t_vec, anomaly_score(worst_m,:), 'r-','LineWidth',2.5,...
    'DisplayName',sprintf('Worst m=%d',worst_m));
plot(t_vec, anomaly_score(best_m,:),  'b-','LineWidth',2.5,...
    'DisplayName',sprintf('Best m=%d',best_m));
plot(t_vec, mean(anomaly_score,1),    'k-','LineWidth',2,...
    'DisplayName','Ensemble mean');
hold off;
xlabel('Forecast step $k$','Interpreter','latex','FontSize',12);
ylabel('$\|\delta\kappa_k^i\|_\mathcal{H}^2$',...
    'Interpreter','latex','FontSize',12);
title('RKHS anomaly score (worst/best)',...
    'FontSize',11);
legend('Location','best','FontSize',9,'Interpreter','latex');
grid on;

sgtitle(['Dual amplitudes $\delta\alpha_j^i = '...
    '\langle\varphi_j,\kappa_{x_0^i}-\bar\kappa_{x_0}\rangle_\mathcal{H}$'],...
    'Interpreter','latex','FontSize',14,'FontWeight','bold');



%% Plot 4 -1 :===========================
% Alpha_pairwise: M × M × N_modes  — full pairwise
% for M=100, N_modes=50: 100×100×50 = 500,000 numbers — fine
k_show = 10;

Alpha_pairwise = zeros(M, M, N_modes);
for i = 1:M
    for l = 1:M
        Alpha_pairwise(i,l,:) = Phi0(i,:) - Phi0(l,:);
    end
end

% Vectorized — much cleaner
% Alpha_pairwise(i,l,j) = Phi0(i,j) - Phi0(l,j)
% This is a 3D outer difference — no loop needed
Alpha_pairwise = permute(Phi0, [1,3,2]) - permute(Phi0, [3,1,2]);
% M × M × N_modes

% Pairwise anomaly score matrix at each time k
% score_pairwise: M × M × (K_forecast+1)
weights_all = abs(KEs').^2;   % 1 × N_modes

score_pairwise = zeros(M, M, K_forecast+1);
for k_t = 0:K_forecast
    weights = abs(KEs).^(2*k_t);    % N_modes × 1
    % |Alpha_pairwise|^2 * weights: M × M
    score_pairwise(:,:,k_t+1) = ...
        sum(abs(Alpha_pairwise).^2 .* permute(weights,[3,2,1]), 3);
end


%% ============================================================
%% PLOT 4 REVISED: Full pairwise dual amplitudes
%% ============================================================
figure('Color','w','Position',[100 100 1200 500]);

% Panel 1: pairwise score matrix at k=0 (initial)
subplot(1,3,1);
imagesc(1:M, 1:M, log(score_pairwise(:,:,1)));
colorbar; colormap(jet); axis square;
xlabel('Member $\ell$','Interpreter','latex','FontSize',12);
ylabel('Member $i$',   'Interpreter','latex','FontSize',12);
title('$\|\delta\kappa_0^{(i,\ell)}\|^2_\mathcal{H}$  at $k=0$',...
    'Interpreter','latex','FontSize',12);

% Panel 2: pairwise score matrix at k=k_show
subplot(1,3,2);
imagesc(1:M, 1:M, log(score_pairwise(:,:,k_show+1)));
colorbar; colormap(jet); axis square;
xlabel('Member $\ell$','Interpreter','latex','FontSize',12);
ylabel('Member $i$',   'Interpreter','latex','FontSize',12);
title(sprintf('$\\|\\delta\\kappa_{%d}^{(i,\\ell)}\\|^2_\\mathcal{H}$  at $k=%d$',...
    k_show, k_show),'Interpreter','latex','FontSize',12);

% Panel 3: row-sum = total anomaly of member i against all others
subplot(1,3,3);
total_score_k0   = sum(score_pairwise(:,:,1),        2);  % M × 1
total_score_kend = sum(score_pairwise(:,:,end),       2);  % M × 1

[~, worst_m] = max(total_score_kend);
[~, best_m]  = min(total_score_kend);

hold on;
bar(1:M, total_score_kend, 'FaceColor',[0.7 0.7 0.9], 'EdgeColor','none');
bar(worst_m, total_score_kend(worst_m), 'FaceColor','r');
bar(best_m,  total_score_kend(best_m),  'FaceColor','b');
hold off;
xlabel('Ensemble member $i$','Interpreter','latex','FontSize',12);
ylabel('$\sum_\ell \|\delta\kappa_k^{(i,\ell)}\|^2_\mathcal{H}$',...
    'Interpreter','latex','FontSize',11);
%title('Total pairwise anomaly per member',...
   % 'FontSize',11);
grid on;

sgtitle(['Pairwise dual amplitudes $\delta\alpha_j^{(i,\ell)} = '...
    '\langle\varphi_j,\,\kappa_{x_0^i}-\kappa_{x_0^\ell}'...
    '\rangle_{\mathcal{H}}$'],...
    'Interpreter','latex','FontSize',13,'FontWeight','bold'); 

fprintf('Worst member (max total pairwise score): m=%d\n', worst_m);
fprintf('Best  member (min total pairwise score): m=%d\n', best_m);
box on;




