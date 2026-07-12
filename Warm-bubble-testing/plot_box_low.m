% plot_box_low-order_surrogate_model.m
%% --- 1. Setup Variables and Matrices ---
KEs_sorted  = KEs(sort_idx);
KEFs_sorted = KEFs(:, sort_idx);
%KMs_sorted  = KMs(:, sort_idx);

r_values = [5, 15, 20, 30, 50, 70, 100, 120, 150]; 
num_r = length(r_values);
total_snapshots = size(Ya_all, 1); 

% Two matrices: one for empirical residuals, one for theoretical bounds
log_residual_matrix = zeros(total_snapshots, num_r);
log_bound_matrix    = zeros(total_snapshots, num_r);

% Global baseline parameter for the theorem's normalization
total_eigen_mass = sum(abs(KEs_sorted));

%% 2. Unified One-Shot Matrix Loop
for r_idx = 1:num_r
    r = r_values(r_idx);
    
    % Slicing variables cleanly based on rank r
    KEs_r   = KEs_sorted(1:r);       
    %KMs_r   = KMs_sorted(:, 1:r);
    KEFs_r  = KEFs_sorted(:, 1:r);

    % 2. THE CRITICAL MATHEMATICAL FIX:
    % and pseudo-inverse relationships match perfectly!
     KMs_r = Xa_all * pinv(KEFs_r.');  %KMs_r = Xa_all / (KEFs_r.');
  
    % --- 1. EMPIRICAL RESIDUAL CALCULATION ---
    Ya_pred_r = real(KMs_r * diag(KEs_r) * KEFs_r.'); 
    
    true_norms = sqrt(sum(Ya_all.^2, 2));
    error_norms = sqrt(sum((Ya_all - Ya_pred_r).^2, 2));
    res_norm = error_norms ./ true_norms;
    res_norm(res_norm < 1e-12) = 1e-12; % Protect against log10(0)
    
    log_residual_matrix(:, r_idx) = log10(res_norm);
    
    % --- 2. THEORETICAL ERROR BOUND DISTRIBUTION CALCULATION ---
    % Discarded tail energy from your theorem bounds
    unmodeled_eigen_mass = sum(abs(KEs_sorted(r+1:end)));
    
    % Snapshot-specific scaling variation (if your theorem scales by snapshot norm)
    % Here it projects the bound dynamically across your snapshot spectrum
    for t = 1:total_snapshots
        % Fine bound equation derived from your dual-space theorem
        snap_bound = unmodeled_eigen_mass / total_eigen_mass;
        if snap_bound < 1e-12, snap_bound = 1e-12; end
        log_bound_matrix(t, r_idx) = log10(snap_bound);
    end
end

%% --- 3. Plotting the Combined Box-Violin Visual ---
figure('Position', [100, 100, 800, 500]);
hold on;

% Color configurations for your paper layout
residual_color = [0.2 0.4 0.8]; % Blue for Residuals
bound_color    = [0.8 0.2 0.2]; % Red/Orange for Bounds

% 1. Plot Theoretical Error Bounds as Violin Distributions
if  exist('violinplot', 'file')
    % If using MATLAB's native violinplot (introduced in recent versions)
    v = violinplot(log_bound_matrix, string(r_values), 'ViolinColor', bound_color, 'EdgeColor', 'none');
    for k = 1:length(v), v(k).BoxPlot.Visible = 'off'; v(k).WhiskerPlot.Visible = 'off'; end
else
    % Fallback: Plot bounds as a subtle shaded distribution curve using standard patch
    for r_idx = 1:num_r
        [f, xi] = ksdensity(log_bound_matrix(:, r_idx));
        f = f / max(f) * 0.3; % Scale width to fit nicely between x-ticks
        patch(r_idx + f, xi, bound_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        patch(r_idx - f, xi, bound_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
end

% 2. Overlay Empirical Residuals as crisp Boxplots shifted slightly for clarity
boxplot(log_residual_matrix, 'Labels', cellstr(string(r_values)), 'Widths', 0.3, 'Colors', residual_color);

%% --- 4. Professional Formatting & Legend ---
grid on;
set(gca, 'GridLineStyle', ':', 'LineWidth', 1.2, 'FontSize', 12);
xlabel('Surrogate Model Order (Rank $r$)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Log Scale $\log_{10}(\text{Value})$', 'Interpreter', 'latex', 'FontSize', 14);
title('Comparison: Empirical Residuals (Box) vs. Theorem Fine Error Bounds (Violin)', ...
    'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');

% Create custom legend handles to keep it clean
h1 = plot(NaN, NaN, 's', 'Color', residual_color, 'MarkerFaceColor', residual_color, 'LineWidth', 2);
h2 = fill(NaN, NaN, bound_color, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
legend([h1, h2], {'Empirical Snapshot Residuals', 'Theorem Fine Error Bounds'}, ...
    'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 11);

hold off;






%% ====================
%% --- 1. Compute the Specific Low-Rank Slices for Comparison ---

%  Ya_pred = KMs * diag(KEs) * KEFs.';
% Extract rank r = 15 components
r1 = 120;
KEs_r1   = KEs_sorted(1:r1);
KEFs_r1  = KEFs_sorted(:, 1:r1);

% Re-compute KMs_r1 algebraically from the sliced subspace
 KMs_r1 = Xa_all * pinv(KEFs_r1.'); %KMs_r1   = Xa_all / (KEFs_r1.');

Ya_pred_r1_matrix = real(KMs_r1 * diag(KEs_r1) * KEFs_r1.');

% Extract rank r = 50 components
r2 = 70;

KEs_r2   = KEs_sorted(1:r2);
KEFs_r2  = KEFs_sorted(:, 1:r2);

% Re-compute KMs_r1 algebraically from the sliced subspace
 KMs_r2 = Xa_all * pinv(KEFs_r2.');   %% KMs_r2   = Xa_all / (KEFs_r2.');

Ya_pred_r2_matrix = real(KMs_r2 * diag(KEs_r2) * KEFs_r2.');

%% --- 2. Combined Trajectory Comparison Plot ---
figure('Position', [100, 100, 850, 480]); 
hold on;

% 1. Create the shaded background regions FIRST so the lines sit on top
colors = num2cell(lines(M), 2); % Generate M distinct colors for shading
for i = 1:M
    % Calculate the start and end point for this ensemble
    start_pt = (i-1)*T_snap + 1;
    end_pt = i*T_snap;
    
    % Draw a shaded background block for the current ensemble
    xregion(start_pt, end_pt, 'FaceColor', colors{i}, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    
    % Add a text label near the top of each shaded block
    text(start_pt + T_snap/2, max(real(Ya_all(indx_spec,:))) * 0.9, ...
        sprintf('Ens# %d', i), 'HorizontalAlignment', 'center', ...
        'FontSize', 20, 'FontWeight', 'bold', 'Color', colors{i}*0.7);
end

% 2. Plot the baseline trajectories (True vs Full Prediction)
p1 = plot(Ya_all(indx_spec,:), 'b.-', 'LineWidth', 1.5, 'DisplayName', 'True State');
p2 = plot(real(Ya_pred(indx_spec,:)), 'r-', 'LineWidth', 1.2, 'DisplayName', 'Full Model');

% 3. Overlay the Low-Rank Surrogates for comparison
% Using dashed lines and distinctive colors to keep the plot readable
p3 = plot(Ya_pred_r1_matrix(indx_spec,:), 'g--', 'LineWidth', 1.1, ...
    'DisplayName', "Surrogate (r = " + r1 + ")");

p4 = plot(Ya_pred_r2_matrix(indx_spec,:), 'm--', 'LineWidth', 1.1, ...
    'DisplayName', "Surrogate (r = " + r2 + ")");

% 4. Clean up the plot aesthetics
grid on;
set(gca, 'GridLineStyle', ':', 'LineWidth', 1.1);
xlim([1, size(Ya_all, 2)]);
xlabel('Snapshot (All Ensemble Data Points)', 'FontSize', 12);
ylabel(sprintf('Some Index of Specific Location %d', indx_spec), 'FontSize', 12);
title(sprintf('Ensemble Trajectory Comparisons (State Variable %d)', indx_spec), 'FontSize', 13, 'FontWeight', 'bold');

% Include all lines cleanly in the legend
legend([p1, p2, p3, p4], 'Location', 'best', 'FontSize', 10);
ax = gca; ax.FontSize = 12; box on;
hold off;




%% rechek the abnormal ?? 


%% =====================================================================
%% 1. FIXED LOAD AND INITIAL PARAMETERS
%% =====================================================================
Xa_all = EnVfull_Xa_36;  % Matrix size: [3880 x 360]
Ya_all = EnVfull_Ya_36;  % Matrix size: [3880 x 360]

grid_points = size(Ya_all, 1);           % 3880 spatial states
T_snap = 36;                             % Snapshots per ensemble
M = 10;                                  % Total number of ensembles

% Safe Rank Limit for Dual KMD Subspace Tracking
max_rank = size(KEFs, 2);                % 360
r = min(100, max_rank);                  % Chosen truncation rank

% --- ONLY SORT KEs AND KEFs ---
KEs_sorted  = KEs(sort_idx);
KEFs_sorted = KEFs(:, sort_idx);

% --- CORRECT ALGEBRAIC SUBSPACE KMs RECONSTRUCTION ---
KEs_r  = KEs_sorted(1:r);
KEFs_r = KEFs_sorted(:, 1:r);
KMs_r  = Xa_all / (KEFs_r.');            % Recompute fresh to prevent scrambled columns

%% =====================================================================
%% 2. METHOD 1: COMPUTE RAW PHYSICAL ENSEMBLE MEAN & DEVIATIONS
%% =====================================================================
Xa_EnMean = zeros(grid_points, T_snap);
Ya_EnMean = zeros(grid_points, T_snap);

for i = 1:M
    col_indices = (i-1)*T_snap + (1:T_snap);
    Xa_EnMean = Xa_EnMean + Xa_all(:, col_indices);
    Ya_EnMean = Ya_EnMean + Ya_all(:, col_indices);
end
Xa_EnMean = Xa_EnMean / M;
Ya_EnMean = Ya_EnMean / M;

raw_deviation_tracks = zeros(T_snap, M);
for i = 1:M
    col_indices = (i-1)*T_snap + (1:T_snap);
    for t = 1:T_snap
        raw_deviation_tracks(t, i) = norm(Ya_all(:, col_indices(t)) - Ya_EnMean(:, t), 2);
    end
end
cumulative_raw_scores = sqrt(sum(raw_deviation_tracks.^2, 1))';
[~, raw_abnormal_idx] = max(cumulative_raw_scores);
[~, raw_best_idx]     = min(cumulative_raw_scores);

%% =====================================================================
%% 3. METHOD 2: ANALYTICAL DUAL KMD SUBSPACE TRAJECTORY TRACKING
%% =====================================================================
initial_indices = (0:(M-1)) * T_snap + 1;
KEFs_init = KEFs_r(initial_indices, :);  % Size: [M x r]
mean_c0 = mean(KEFs_init, 1);            % Size: [1 x r]

kmd_subspace_deviation_norms = zeros(T_snap, M);
for t = 1:T_snap
    lambda_term = (KEs_r.').^(t-1);      % Analytical time propagation factor
    for i = 1:M
        c0_deviation = KEFs_init(i, :) - mean_c0;
        propagated_dev_features = c0_deviation .* lambda_term;
        
        % Map back to physical space via the dynamically correct KMs_r matrix
        d_t_i = real(KMs_r * propagated_dev_features.');
        kmd_subspace_deviation_norms(t, i) = norm(d_t_i, 2);
    end
end
cumulative_kmd_scores = sqrt(sum(kmd_subspace_deviation_norms.^2, 1))';
[~, kmd_abnormal_idx] = max(cumulative_kmd_scores);
[~, kmd_best_idx]     = min(cumulative_kmd_scores);

%% =====================================================================
%% 4. NEW PLOT: SIDE-BY-SIDE TRAJECTORY TRACK COMPARISON
%% =====================================================================
figure('Name', 'Trajectory Deviation Comparison', 'Position', [100, 100, 1100, 420]);

% --- LEFT SUBPLOT: RAW DATA DEVIATIONS ---
subplot(1, 2, 1); hold on; grid on;
for i = 1:M
    if i ~= raw_abnormal_idx && i ~= raw_best_idx
        plot(1:T_snap, raw_deviation_tracks(:, i), 'Color', [0.8 0.8 0.8], 'LineWidth', 1.2);
    end
end
p1 = plot(1:T_snap, raw_deviation_tracks(:, raw_best_idx), 'Color', [0.2 0.4 0.8], 'LineWidth', 2);
p2 = plot(1:T_snap, raw_deviation_tracks(:, raw_abnormal_idx), 'Color', [0.8 0.2 0.2], 'LineWidth', 3);
set(gca, 'GridLineStyle', ':', 'FontSize', 11);
xlim([1, T_snap]);
xlabel('Ensemble Time Step (t)', 'FontSize', 12);
ylabel('Physical Deviation Norm \|y_t^i - y_{mean,t}\|', 'FontSize', 12);
title('Method 1: Raw Physical Space Tracks', 'FontSize', 12, 'FontWeight', 'bold');
legend([p1, p2], {['Best (Idx ', num2str(raw_best_idx), ')'], ['Abnormal (Idx ', num2str(raw_abnormal_idx), ')']}, 'Location', 'best');

% --- RIGHT SUBPLOT: DUAL KMD SUBSPACE DEVIATIONS ---
subplot(1, 2, 2); hold on; grid on;
for i = 1:M
    if i ~= kmd_abnormal_idx && i ~= kmd_best_idx
        plot(1:T_snap, kmd_subspace_deviation_norms(:, i), 'Color', [0.8 0.8 0.8], 'LineWidth', 1.2);
    end
end
p3 = plot(1:T_snap, kmd_subspace_deviation_norms(:, kmd_best_idx), 'Color', [0.2 0.4 0.8], 'LineWidth', 2);
p4 = plot(1:T_snap, kmd_subspace_deviation_norms(:, kmd_abnormal_idx), 'Color', [0.8 0.2 0.2], 'LineWidth', 3);
set(gca, 'GridLineStyle', ':', 'FontSize', 11);
xlim([1, T_snap]);
xlabel('Ensemble Time Step (t)', 'FontSize', 12);
ylabel('Subspace Deviation Norm \|d_t^i\|', 'FontSize', 12);
title(sprintf('Method 2: Analytical Dual KMD Subspace (r=%d)', r), 'FontSize', 12, 'FontWeight', 'bold');
legend([p3, p4], {['Best (Idx ', num2str(kmd_best_idx), ')'], ['Abnormal (Idx ', num2str(kmd_abnormal_idx), ')']}, 'Location', 'best');

%% =====================================================================
%% 5. NEW PLOT: SIDE-BY-SIDE CUMULATIVE ENERGY BAR CHARTS
%% =====================================================================
figure('Name', 'Cumulative Deviation Comparison', 'Position', [150, 150, 1000, 400]);

% --- LEFT BAR CHART: RAW ENERGY ---
subplot(1, 2, 1);
bar(1:M, cumulative_raw_scores, 'FaceColor', [0.65 0.65 0.65], 'EdgeColor', 'none'); hold on;
bar(raw_best_idx, cumulative_raw_scores(raw_best_idx), 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none');
bar(raw_abnormal_idx, cumulative_raw_scores(raw_abnormal_idx), 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none');
grid on; set(gca, 'GridLineStyle', ':', 'XTick', 1:M, 'FontSize', 11);
xlabel('Ensemble Member Index i', 'FontSize', 12);
ylabel('Integrated Raw Physical Energy', 'FontSize', 12);
title('Raw Space Integrated Total', 'FontSize', 12, 'FontWeight', 'bold');

% --- RIGHT BAR CHART: SUBSPACE ENERGY ---
subplot(1, 2, 2);
bar(1:M, cumulative_kmd_scores, 'FaceColor', [0.65 0.65 0.65], 'EdgeColor', 'none'); hold on;
bar(kmd_best_idx, cumulative_kmd_scores(kmd_best_idx), 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none');
bar(kmd_abnormal_idx, cumulative_kmd_scores(kmd_abnormal_idx), 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none');
grid on; set(gca, 'GridLineStyle', ':', 'XTick', 1:M, 'FontSize', 11);
xlabel('Ensemble Member Index i', 'FontSize', 12);
ylabel('Integrated Subspace Energy', 'FontSize', 12);
title('Subspace Projected Integrated Total', 'FontSize', 12, 'FontWeight', 'bold');




%%%

%% 1. LOAD AND SETUP PARAMETERS
% Xa_all and Ya_all are size [3880 x 360] (3880 grid points, 360 total snapshots)
Xa_all = EnVfull_Xa_36;  
Ya_all = EnVfull_Ya_36;

T_snap = 36;                             % Snapshots per ensemble
M = 10;                                  % Total number of ensembles
grid_points = size(Ya_all, 1);           % 3880 spatial states

%% 2. COMPUTE THE RAW ENSEMBLE MEAN TRAJECTORY SIGNALS
% We isolate each ensemble block of 36 columns and average them together 
% to get one single reference mean trajectory matrix of size [3880 x 360]
Xa_EnMean = zeros(grid_points, T_snap);
Ya_EnMean = zeros(grid_points, T_snap);

for i = 1:M
    col_indices = (i-1)*T_snap + (1:T_snap);
    Xa_EnMean = Xa_EnMean + Xa_all(:, col_indices);
    Ya_EnMean = Ya_EnMean + Ya_all(:, col_indices);
end
Xa_EnMean = Xa_EnMean / M;
Ya_EnMean = Ya_EnMean / M;

%% 3. TRACK THE RAW PHYSICAL TRAJECTORY DEVIATIONS
% Calculate the norm deviation of each ensemble from the raw mean signal at every time step
raw_deviation_tracks = zeros(T_snap, M);

for i = 1:M
    col_indices = (i-1)*T_snap + (1:T_snap);
    
    for t = 1:T_snap
        % Current snapshot column for ensemble i at time step t
        current_snapshot = Ya_all(:, col_indices(t));
        
        % Mean snapshot column at time step t
        mean_snapshot = Ya_EnMean(:, t);
        
        % True physical Euclidean distance over all 3880 states
        raw_deviation_tracks(t, i) = norm(current_snapshot - mean_snapshot, 2);
    end
end

% Compute the cumulative space-time energy deviation score for each ensemble
cumulative_deviation_scores = sqrt(sum(raw_deviation_tracks.^2, 1))';

% Automatically identify the best and worst ensemble indexes based on your logic
[~, abnormal_idx] = max(cumulative_deviation_scores);
[~, best_idx]     = min(cumulative_deviation_scores);

%% 4. PLOT THE ACTUAL PHYSICAL TRAJECTORY TRACKS
figure('Name', 'Raw Physical Trajectory Deviations', 'Position', [100, 100, 620, 420]);
hold on; grid on;

% Plot normal cluster tracks in light grey
for i = 1:M
    if i ~= abnormal_idx && i ~= best_idx
        plot(1:T_snap, raw_deviation_tracks(:, i), 'Color', [0.75 0.75 0.75], 'LineWidth', 1.2);
    end
end

% Highlight the best nominal track (closest to the raw mean signal) in blue
p_best = plot(1:T_snap, raw_deviation_tracks(:, best_idx), 'Color', [0.2 0.4 0.8], 'LineWidth', 2);

% Highlight the abnormal outlier track (largest deviation from the raw mean signal) in red
p_abnormal = plot(1:T_snap, raw_deviation_tracks(:, abnormal_idx), 'Color', [0.8 0.2 0.2], 'LineWidth', 3);

set(gca, 'GridLineStyle', ':', 'FontSize', 12);
xlim([1, T_snap]);
xlabel('Ensemble Time Step (t)', 'FontSize', 13);
ylabel('Physical Deviation Norm ||y_t^i - y_{mean,t}||', 'FontSize', 13);
title('Physical Trajectory Deviation Tracks from Raw Ensemble Mean', 'FontSize', 13, 'FontWeight', 'bold');

legend([p_best, p_abnormal], ...
    {['Best Nominal Track (Index ', num2str(best_idx), ')'], ...
     ['Abnormal Outlier Track (Index ', num2str(abnormal_idx), ')']}, ...
    'Location', 'best');
hold off;

%% 5. PLOT THE CUMULATIVE DEVIATION BAR CHART
figure('Name', 'Cumulative Raw Deviations', 'Position', [740, 100, 520, 420]);
b = bar(1:M, cumulative_deviation_scores, 'FaceColor', [0.65 0.65 0.65], 'EdgeColor', 'none');
hold on;

% Color highlights for index indicators matching the tracks
bar(best_idx, cumulative_deviation_scores(best_idx), 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none');
bar(abnormal_idx, cumulative_deviation_scores(abnormal_idx), 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none');

grid on;
set(gca, 'GridLineStyle', ':', 'XTick', 1:M, 'FontSize', 12);
xlabel('Ensemble Member Index i', 'FontSize', 13);
ylabel('Integrated Physical Deviation Energy', 'FontSize', 13);
title('Total Cumulative Space-Time Deviation Summary', 'FontSize', 13, 'FontWeight', 'bold');
hold off;