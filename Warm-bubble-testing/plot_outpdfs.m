%%% plot_pdf_outs.m
%% =============Pdf-output
figure('Position', [100, 100, 1300, 1100]);
numRows = 3;
numCols = 3;

% Use 'compact' spacing to keep tiles close while accommodating individual colorbars
ttt = tiledlayout(numRows, numCols, 'Padding', 'compact', 'TileSpacing', 'compact');

% Skip the complex conjugates to only plot unique spatial patterns!
modesToPlot = 1:2:17; 

for k = 1:length(modesToPlot)
    % This is the TRUE sorted index (1, 3, 5, 7, ..., 17)
    mode_index_in_sorted = modesToPlot(k);
    
    % Extract the spatial mode column vector
    KMode_i = abs(KMs(:, mode_index_in_sorted)); 
    
    % Reshape perfectly back to your physical grid layout
    mode_i_reshaped = reshape(KMode_i, grid_size); 
    
    % Automatically jumps to the next available slot in the 3x3 layout
    ax = nexttile;
    
    imagesc(data.y, data.z, mode_i_reshaped); 
    colormap(ax, brighten(redblueTecplot(21), -0.55));
    
    % Add individual colorbar for each subplot
    cb = colorbar('Location', 'eastoutside');
    cb.LineWidth = 0.75;
    
    % Maintain physical coordinates tracking orientation (z going upwards)
    axis xy;
    axis equal; 
    xlim([0 2e4]);
    ylim([0 2e4]);
    
    set(gca, 'FontSize', 13);    
    xlabel("y", 'FontSize', 12); 
    ylabel("z", 'FontSize', 12);
    % Added 'Interpreter', 'latex' here so y and z render in LaTeX style
    %xlabel("$y$", 'Interpreter', 'latex', 'FontSize', 10); 
    %ylabel("$z$", 'Interpreter', 'latex', 'FontSize', 10); 
    
    % 1. Extract the current complex eigenvalue
    lambda_j = KEs(mode_index_in_sorted);
    real_part = real(lambda_j);
    imag_part = imag(lambda_j);
    
    % 2. Clean LaTeX String (sprintf ONLY formats the text string)
   % title_str = sprintf('$\\textrm{Mode }%d: \\lambda_{%d} = %.2f %+.2f i$', ...
   %     mode_index_in_sorted, mode_index_in_sorted, real_part, imag_part);

   % 2. Clean LaTeX String (Removed "Mode %d:" prefix)
    title_str = sprintf('$\\lambda_{%d} = %.2f %+.2f i$', ...
        mode_index_in_sorted, real_part, imag_part);
    
    % 3. Plot with LaTeX Interpreter and FontSize properly passed to title()
    title(title_str, 'Interpreter', 'latex', 'FontSize', 11);
    
    %% --- 4. Direct Coordinate Marker Overlay ---
    highlight_y = data.y(highlight_col); 
    highlight_z = data.z(highlight_row); 
    
    hold on;
    plot(highlight_y, highlight_z, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow');
    hold off;
end

set(gcf, 'Renderer', 'painters');

% Export the complete 3x3 layout as a single clean vector PDF
exportgraphics(gcf, 'ensemble_modes_3x3_individual_cbar.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none');



%% =============================================
%% 4. PLOT THE ACTUAL PHYSICAL TRAJECTORY TRACKS
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

fig2 = figure('Name', 'Raw Physical Trajectory Deviations', 'Position', [100, 100, 620, 420]);
hold on; grid on;

% Plot normal cluster tracks in light grey
for i = 1:M
    if i ~= abnormal_idx && i ~= best_idx
        plot(1:T_snap, raw_deviation_tracks(:, i), 'Color', [0.75 0.75 0.75], 'LineWidth', 1.2);
    end
end

% Highlight the best nominal track (closest to the raw mean signal) in blue
p_best = plot(1:T_snap, raw_deviation_tracks(:, best_idx), 'Color', [0.2 0.4 0.8], 'LineWidth', 3);

% Highlight the abnormal outlier track (largest deviation from the raw mean signal) in red
p_abnormal = plot(1:T_snap, raw_deviation_tracks(:, abnormal_idx), 'Color', [0.8 0.2 0.2], 'LineWidth', 3);

% Axis formatting
ax = gca;
ax.GridLineStyle = ':';
ax.FontSize = 20;                     % x and y tick label font size
ax.TickLabelInterpreter = 'latex';    % LaTeX-render the tick numbers too
xlim([1, T_snap]);
box on
set(gca, 'LineWidth', 1.5)

xlabel('Snapshot $t$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Physical Deviation Norm $\|y_t^i - \bar{y}_t\|$', 'Interpreter', 'latex', 'FontSize', 20);
title('Trajectory Deviation', ...
    'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold');

legend([p_best, p_abnormal], ...
    {['Best Nominal (Index ', num2str(best_idx), ')'], ...
     ['Abnormal Track (Index ', num2str(abnormal_idx), ')']}, ...
    'Interpreter', 'latex', 'Location', 'best');

hold off;
set(gcf, 'Renderer', 'painters');

exportgraphics(fig2, 'abnomral_detetion_physical.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none');

%%    out-of samples pdf 

% ==========================================
% ---- FIGURE 1: In-Sample (Training Fit) ----

r_values = [15, 30, 45, 60, 75, 90, 105, 120, 150];
num_r    = length(r_values);

fig1 = figure('Position', [100, 100, 600, 450]);

if use_violin
    violinplot(log_residual_in, string(r_values), 'ViolinColor', [0.20, 0.40, 0.75]);
else
    h = boxplot(log_residual_in, 'Colors', [0.20, 0.40, 0.75], 'Symbol', 'o');
    set(h, 'LineWidth', 1.75);

    outliers = findobj(h, 'Tag', 'Outliers');
    set(outliers, 'MarkerSize', 6, 'LineWidth', 1.5);

    n_cols_in = size(log_residual_in, 2);
    if n_cols_in ~= num_r
        warning(['Mismatch: log_residual_in has ', num2str(n_cols_in), ...
                 ' columns but r_values has ', num2str(num_r), ' entries.']);
    end

    set(gca, 'XTick', 1:n_cols_in);
    if n_cols_in <= num_r
        set(gca, 'XTickLabel', string(r_values(1:n_cols_in)));
    else
        labels = [string(r_values), repmat("?", 1, n_cols_in - num_r)];
        set(gca, 'XTickLabel', labels);
    end
end

grid on;

% --- Axis tick label font size (x and y) ---
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);

% --- Axis label / title font size ---
xlabel('No. of Modes', 'FontWeight', 'FontSize', 20);
ylabel('Relative Error', 'FontWeight', 'FontSize', 20);
title('Error bound', 'FontWeight', 'bold', 'FontSize', 20);

set(gcf, 'Renderer', 'painters');

exportgraphics(fig1, 'error_bounds_in_samples.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none');




% ==========================================
% ---- FIGURE 2: Out-of-Sample -----------
% ==========================================
fig2 = figure('Position', [750, 100, 600, 450]);

if n_test_snap > 1
    if use_violin
        violinplot(log_residual_out, string(r_values), 'ViolinColor', [0.85, 0.30, 0.25]);
    else
        boxplot(log_residual_out, 'Colors', [0.85, 0.30, 0.25], 'Symbol', 'o');
        
        n_cols_out = size(log_residual_out, 2);
        set(gca, 'XTick', 1:n_cols_out);
        if length(r_values) >= n_cols_out
            set(gca, 'XTickLabel', r_values(1:n_cols_out));
        else
            set(gca, 'XTickLabel', r_values);
        end
    end
else
    % Only one out-of-sample target -> plot as a line
    plot(r_values, log_residual_out, '-o', 'Color', [0.85, 0.30, 0.25], ...
        'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', [0.85, 0.30, 0.25]);
    grid on;
end

set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
xlabel('Retained Modes (r)', 'FontWeight', 'bold');
ylabel('log_{10}(Relative Error)', 'FontWeight', 'bold');
title('Out-of-Sample (Ensemble Mean, Xb/Yb)', 'FontWeight', 'bold');

set(gcf, 'Renderer', 'painters');
exportgraphics(fig2, 'out_of_sample_error.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none');



%% -------------
%  Ya_pred = KMs * diag(KEs) * KEFs.';
% Extract rank r = 15 components
% r1 = 120;
% KEs_r1   = KEs_sorted(1:r1);
% KEFs_r1  = KEFs_sorted(:, 1:r1);
% 
% % Re-compute KMs_r1 algebraically from the sliced subspace
%  KMs_r1 = Xa_all * pinv(KEFs_r1.'); %KMs_r1   = Xa_all / (KEFs_r1.');
% 
% Ya_pred_r1_matrix = real(KMs_r1 * diag(KEs_r1) * KEFs_r1.');
% 
% % Extract rank r = 50 components
% r2 = 70;
% 
% KEs_r2   = KEs_sorted(1:r2);
% KEFs_r2  = KEFs_sorted(:, 1:r2);
% 
% % Re-compute KMs_r1 algebraically from the sliced subspace
%  KMs_r2 = Xa_all * pinv(KEFs_r2.');   %% KMs_r2   = Xa_all / (KEFs_r2.');
% 
% Ya_pred_r2_matrix = real(KMs_r2 * diag(KEs_r2) * KEFs_r2.');

%% --- 2. Combined Trajectory Comparison Plot ---
figure('Position', [100, 100, 850, 480]);
hold on;

% 1. Create the shaded background regions FIRST so the lines sit on top
colors = num2cell(lines(M), 2);
for i = 1:M
    start_pt = (i-1)*T_snap + 1;
    end_pt = i*T_snap;
    xregion(start_pt, end_pt, 'FaceColor', colors{i}, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    text(start_pt + T_snap/2, max(real(Ya_all(indx_spec,:))), ...   % <-- moved up: 0.9 -> 0.97
        sprintf('Ens# %d', i), 'HorizontalAlignment', 'center', ...
        'FontSize', 15, 'FontWeight', 'bold', 'Color', colors{i}*0.7);
end

% --- Strong, saturated, well-separated colors ---
c1 = [0.00, 0.00, 0.85];
c2 = [0.85, 0.00, 0.00];
c3 = [0.00, 0.60, 0.00];
c4 = [0.75, 0.00, 0.75];

n_pts = size(Ya_all, 2);
mark_step = max(round(n_pts / 40), 1);

p2 = plot(real(Ya_pred(indx_spec,:)), '-.', 'Color', c2, 'LineWidth', 1.2, ...
    'Marker', 's', 'MarkerSize', 9, 'MarkerFaceColor', c2, 'MarkerEdgeColor', 'k', ...
    'MarkerIndices', 1:mark_step:n_pts, 'DisplayName', 'Full Model');

p3 = plot(Ya_pred_r1_matrix(indx_spec,:), '--', 'Color', c3, 'LineWidth', 2.0, ...
    'Marker', '^', 'MarkerSize', 10, 'MarkerFaceColor', c3, 'MarkerEdgeColor', 'k', ...
    'MarkerIndices', 1:mark_step:n_pts, 'DisplayName', "Surrogate (r = " + r1 + ")");

p4 = plot(Ya_pred_r2_matrix(indx_spec,:), ':', 'Color', c4, 'LineWidth', 2.0, ...
    'Marker', 'd', 'MarkerSize', 10, 'MarkerFaceColor', c4, 'MarkerEdgeColor', 'k', ...
    'MarkerIndices', 1:mark_step:n_pts, 'DisplayName', "Surrogate (r = " + r2 + ")");

p1 = plot(Ya_all(indx_spec,:), '-.', 'Color', c1, 'LineWidth', 1.8, ...
    'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', c1, 'MarkerEdgeColor', 'w', ...
    'MarkerIndices', 1:mark_step:n_pts, 'DisplayName', 'True State');

grid on;
set(gca, 'GridLineStyle', ':', 'LineWidth', 1.1);
xlim([1, size(Ya_all, 2)]);
xticks(0:36:360);
xlabel('Snapshot (All Ensemble Data Points)', 'FontSize', 18);
ylabel(sprintf('Specific Location'), 'FontSize', 12);
title(sprintf('Ensemble Trajectory Comparisons'), 'FontSize', 13);

legend([p1, p2, p3, p4], 'Location', 'southeast', 'FontSize', 14);
ax = gca; ax.FontSize = 20; box on;
hold off;
set(gcf, 'Renderer', 'painters');

exportgraphics(gcf, 'ensemble_comparisons_1.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none');