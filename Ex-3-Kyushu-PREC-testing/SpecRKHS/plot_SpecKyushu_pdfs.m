% plot_Spec_kyushu_pdfs.m


%% ================================== Eigenvalues (Plain Plot)
figure('Color', 'w');
hold on;

% 1. Plot Unit Circle
theta = linspace(0, 2*pi, 300);
plot(cos(theta), sin(theta), '--k', 'LineWidth', 1.5);

% 2. Plot Eigenvalues (Re vs Im)
scatter(real(LAM), imag(LAM), 80, [0 0.4470 0.7410], 'o', 'filled', ...
    'LineWidth', 0.8);
%scatter(real(Lambda_res), imag(Lambda_res), 100, [0.8500 0.3250 0.0980], 'x', ...
  %  'LineWidth', 1.5, 'DisplayName', '$\Lambda_{\mathrm{res}}$');

% 3. Axis Formatting & Bounds
axis equal;
axis([-1.5, 1.5, -1.5, 1.5]);
grid on; box on;

% Set axis tick mark font size and LaTeX interpreter
ax = gca;
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';

% Labels, Title & Legend
xlabel('$\mathrm{Re}(\lambda)$', 'Interpreter', 'latex', 'FontSize', 25);
ylabel('$\mathrm{Im}(\lambda)$', 'Interpreter', 'latex', 'FontSize', 25);
title('Eigenvalues', 'Interpreter', 'latex', 'FontSize', 22);
%legend('Interpreter', 'latex', 'FontSize', 16, 'Location', 'best');

% 4. Vector Export
exportgraphics(gcf, 'KE_evals_woRes.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');


%% Additionla-ZOOM-in----lots of lines

figure('Color', 'w', 'Position', [100, 100, 900, 650]);
hold on;

% Define colors
c_dmd  = [0.4, 0.4, 0.4];
c_kedmd  = [0.8500, 0.3250, 0.0980];
c_dual = [0.2, 0.4, 0.8];


% --- 1. Main Plot (All Ensembles + Means) ---
M = size(er1_all, 1);
for m = 1:M
    plot(er2_all(m,:), 'Color', [c_dmd, 0.15],  'LineWidth', 1, 'HandleVisibility', 'off')
    plot(er3_all(m,:), 'Color', [c_kedmd, 0.15],  'LineWidth', 1, 'HandleVisibility', 'off')
    plot(er1_all(m,:), 'Color', [c_dual, 0.15], 'LineWidth', 1, 'HandleVisibility', 'off')
end

plot(mean(er2_all, 1), '--', 'Color', c_dmd,  'LineWidth', 2, 'DisplayName', 'DMD (mean)')
plot(mean(er3_all, 1), '-',  'Color', c_kedmd,  'LineWidth', 2, 'DisplayName', 'KeDMD (mean)')
plot(mean(er1_all, 1), '-.', 'Color', c_dual, 'LineWidth', 2, 'DisplayName', 'DualKoop-RKHS (mean)')

% Main Axes Formatting
set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 25);
grid on; box on;
title('Relative Forecast Errors Comparison', 'Interpreter', 'latex', 'FontSize', 18);
xlabel('Lead Time (Snapshots)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Relative Forecast Error', 'Interpreter', 'latex', 'FontSize', 18);
legend('Interpreter', 'latex', 'FontSize', 16, 'Location', 'northwest');


% % --- 2. Inset Zoom Plot (KeDMD vs DualKoop-RKHS in Log Scale) ---
% % Position: [left, bottom, width, height]
axes('Position', [0.48, 0.22, 0.40, 0.38]);
hold on;

% Plot ensemble lines for KeDMD and DualKoop only
for m = 1:M
    plot(er3_all(m,:), 'Color', [c_kedmd, 0.2],  'LineWidth', 0.8, 'HandleVisibility', 'off');
    plot(er1_all(m,:), 'Color', [c_dual, 0.2], 'LineWidth', 0.8, 'HandleVisibility', 'off');
end

% Plot mean lines
plot(mean(er3_all, 1), '-',  'Color', c_kedmd,  'LineWidth', 2);
plot(mean(er1_all, 1), '-.', 'Color', c_dual, 'LineWidth', 2);

% Inset Formatting (Log Scale)
set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 18);
grid on; box on;
title('\textbf{Zoomed: KeDMD vs DualKoop}', 'Interpreter', 'latex', 'FontSize', 16);

% Auto-tight y-limits around KeDMD and DualKoop data
y_min_zoom = min([er1_all(:); er3_all(:)]);
y_max_zoom = max([er1_all(:); er3_all(:)]);
xlim([1, size(er1_all, 2)]);
ylim([y_min_zoom * 0.9, y_max_zoom * 1.1]);


% --- 3. Export Vector Graphics PDF ---
exportgraphics(gcf, 'ensemble_weather_error_lines.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none');

%% ==============Bar ZOOM-------------


figure('Color', 'w', 'Position', [100, 100, 800, 600]);

c_dmd   = [0.0000 0.4470 0.7410]; % Blue
c_kedmd = [0.8500 0.3250 0.0980]; % Terracotta / Red
c_spec  = [0.9290 0.6940 0.1250]; % Orange


% --- 1. Main Bar Plot ---
b = bar(1:3, rmse_values, 'FaceColor', 'flat', 'BarWidth', 0.5);
b.CData(1,:) = c_dmd;
b.CData(2,:) = c_kedmd;
b.CData(3,:) = c_spec;

% Apply LaTeX interpreter to X-tick labels and axes
set(gca, 'TickLabelInterpreter', 'latex', ...
    'XTick', 1:3, ...
    'XTickLabel', {'\textrm{DMD}', '\textrm{KeDMD}', '\textrm{DualKoop}'});
grid on; box on;
ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 18);
title('Forecast Error Comparison (RMSE)', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 20);


% --- 2. Inset Plot for Small Values (Shifted Slightly Down) ---
% Position: [left, bottom, width, height]
axes('Position', [0.50, 0.32, 0.38, 0.38]); % Changed bottom from 0.45 to 0.32
b_inset = bar(2:3, rmse_values(2:3), 'FaceColor', 'flat', 'BarWidth', 0.5);
b_inset.CData(1,:) = c_kedmd;
b_inset.CData(2,:) = c_spec;

% Apply LaTeX interpreter to inset plot
set(gca, 'TickLabelInterpreter', 'latex', ...
    'XTick', 2:3, ...
    'XTickLabel', {'\textrm{KeDMD}', '\textrm{DualKoop}'});
title('\textbf{Zoomed (KeDMD vs DualKoop)}', 'Interpreter', 'latex', 'FontSize', 11);
grid on; box on;
set(gca, 'FontSize', 14);


% --- 3. Export Vector Graphics PDF ---
exportgraphics(gcf, 'ensemble_overall_rmse_bar_1.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none');



%% ===================Zoom in =============

figure('Color', 'w', 'Position', [100, 100, 900, 650]);
hold on;

% --- 1. Data Preparation ---
N_steps = size(er1_all, 2);
t_vec   = 1:N_steps;
x_poly  = [t_vec, fliplr(t_vec)];

% Compute min-max shaded ranges across ensembles (dim 1)
min1 = min(er1_all, [], 1); max1 = max(er1_all, [], 1);
min2 = min(er2_all, [], 1); max2 = max(er2_all, [], 1);
min3 = min(er3_all, [], 1); max3 = max(er3_all, [], 1);

% Colors
% c_dmd  = [0.0000, 0.4470, 0.7410];
% c_ked  = [0.8500, 0.3250, 0.0980];
% c_dual = [0.9290, 0.6940, 0.1250];

c_dmd   = [0.0000 0.4470 0.7410]; % Blue
c_kedmd = [0.8500 0.3250 0.0980]; % Terracotta / Red
c_spec  = [0.9290 0.6940 0.1250]; % Orange

% --- 2. Main Plot (All 3 Methods) ---
fill(x_poly, [max2, fliplr(min2)], c_dmd,  'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(x_poly, [max3, fliplr(min3)], c_kedmd,  'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(x_poly, [max1, fliplr(min1)], c_dual, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

plot(mean(er2_all, 1), '--', 'Color', c_dmd,  'LineWidth', 2, 'DisplayName', 'DMD (mean)');
plot(mean(er3_all, 1), '-',  'Color', c_kedmd,  'LineWidth', 2, 'DisplayName', 'KeDMD (mean)');
plot(mean(er1_all, 1), '-.', 'Color', c_dual, 'LineWidth', 2, 'DisplayName', 'DualKoop-RKHS (mean)');

% Main Axes Formatting (Main plot x and y tick font size set to 20)
set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 26);
grid on; box on;

title('Relative Forecast Errors Comparison', 'Interpreter', 'latex', 'FontSize', 18);
xlabel('Lead Time (Snapshots)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Relative Forecast Error', 'Interpreter', 'latex', 'FontSize', 18);
legend('Interpreter', 'latex', 'FontSize', 16, 'Location', 'northwest');

% --- 3. Inset Zoom Plot (KeDMD vs DualKoop-RKHS) ---
% Position: [left, bottom, width, height]
axes('Position', [0.50, 0.22, 0.38, 0.38]);
hold on;

% Plot shaded areas for KeDMD and DualKoop
fill(x_poly, [max3, fliplr(min3)], c_kedmd,  'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(x_poly, [max1, fliplr(min1)], c_dual, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Plot mean lines for KeDMD and DualKoop
plot(mean(er3_all, 1), '-',  'Color', c_kedmd,  'LineWidth', 2);
plot(mean(er1_all, 1), '-.', 'Color', c_dual, 'LineWidth', 2);

% Inset Formatting
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
grid on; box on;
title('\textbf{Zoomed: KeDMD vs DualKoop}', 'Interpreter', 'latex', 'FontSize', 12);

% Automatically adjust inset Y-limits based on KeDMD and DualKoop values
y_min_zoom = min([min1, min3]);
y_max_zoom = max([max1, max3]);
xlim([1, N_steps]);
ylim([y_min_zoom * 0.9, y_max_zoom * 1.1]);

% --- 4. Export Vector Graphics PDF ---
exportgraphics(gcf, 'ensemble_shallow_weather_error_1.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none');

%==============Zoom in cdmd ====
%=================================

figure
hold on
p1 = plot(0:steps, log(mean(mean_obs_kedmd_all, 2)), 'linewidth', 3, 'color', [0.8500 0.3250 0.0980]);
p2 = plot(0:steps, log(mean(mean_obs_kmd_all, 2)), '-.','linewidth', 3, 'color', [0.9290 0.6940 0.1250]);
p3 = plot(0:steps, log(mean(mean_obs_exact_all, 2)), '--', 'linewidth', 3, 'color', [0 0.4470 0.7410]);
grid on; box on;

% Set axis tick mark font size and tick label interpreter
ax = gca;
ax.FontSize = 22; % Tick numbers font size
ax.TickLabelInterpreter = 'latex'; % Ensures tick numbers match LaTeX style

title('Mean forecast comparison', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Mean Prediction', 'interpreter', 'latex', 'fontsize', 18)
legend([p3 p1 p2], {'DMD', 'KeDMD', 'DualKoop-RKHS'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'best')

exportgraphics(gcf, 'ensemble_prediction_mean.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')


%=====================
figure
semilogy(0:steps, (mean(Kx0_er_all, 2)), 'linewidth', 3, 'color', [0.9290 0.6940 0.1250])
grid on; box on;

% Set axis tick mark font size and tick label interpreter
ax = gca;
ax.FontSize = 22; % Tick numbers font size
ax.TickLabelInterpreter = 'latex'; % Ensures tick numbers match LaTeX style

title('Relative forecast error for kernel function', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 22)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 22)
legend({'DualKoop-RKHS'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'best')

exportgraphics(gcf, 'ensemble_prediction_kernel_error.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%================
figure
er_obs = mean(abs(obs_exact_all - mean_obs_kmd_all).^2 ./ abs(obs_exact_all).^2, 2);
semilogy(0:steps, (er_obs), 'linewidth', 3, 'color', [0.9290 0.6940 0.1250])
grid on; box on;

% Set axis tick mark font size and tick label interpreter
ax = gca;
ax.FontSize =22; % Tick numbers font size
ax.TickLabelInterpreter = 'latex'; % Ensures tick numbers match LaTeX style

title('Mean relative forecast error', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 22)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 22)
legend({'DualKoop-RKHS'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'best')

exportgraphics(gcf, 'ensemble_error_mean.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')


%%   ===================Revised
figure('Color', 'w', 'Position', [100, 100, 900, 650]);
hold on;

% Define colors
c_dmd  = [0.4, 0.4, 0.4];
c_kedmd  = [0.8500, 0.3250, 0.0980];
c_dual = [0.2, 0.4, 0.8];

% --- 1. Main Plot (All Ensembles + Means) ---
M = size(er1_all, 1);
for m = 1:M
    plot(er2_all(m,:), 'Color', [c_dmd, 0.15],  'LineWidth', 1, 'HandleVisibility', 'off')
    plot(er3_all(m,:), 'Color', [c_kedmd, 0.15],  'LineWidth', 1, 'HandleVisibility', 'off')
    plot(er1_all(m,:), 'Color', [c_dual, 0.15], 'LineWidth', 1, 'HandleVisibility', 'off')
end
plot(mean(er2_all, 1), '--', 'Color', c_dmd,  'LineWidth', 2, 'DisplayName', 'DMD (mean)')
plot(mean(er3_all, 1), '-',  'Color', c_kedmd,  'LineWidth', 2, 'DisplayName', 'KeDMD (mean)')
plot(mean(er1_all, 1), '-.', 'Color', c_dual, 'LineWidth', 2, 'DisplayName', 'DualKoop-RKHS (mean)')

% Main Axes Formatting
set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 25);
grid on; box on;
title('Relative Forecast Errors Comparison', 'Interpreter', 'latex', 'FontSize', 18);
xlabel('Lead Time (Snapshots)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Relative Forecast Error', 'Interpreter', 'latex', 'FontSize', 18);

% Legend moved to bottom-middle, horizontal layout, outside the axes
legend('Interpreter', 'latex', 'FontSize', 16, ...
    'Location', 'southoutside', 'Orientation', 'horizontal');

mainAx = gca;

% --- 2. Inset Zoom Plot (DMD only — deviation-from-mean view), moved UP ---
zoomAx = axes('Position', [0.48, 0.55, 0.40, 0.35]);
hold on;

dmd_mean = mean(er2_all, 1);          % 1 x N
dev = er2_all - dmd_mean;             % M x N, deviation of each member from mean

cmap = gray(M+2);
for m = 1:M
    plot(dev(m,:), 'Color', [cmap(m+1,:), 0.8], 'LineWidth', 1.2, 'HandleVisibility', 'off');
end
plot(zeros(1, size(er2_all,2)), '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2);  % mean = 0 line

set(gca, 'YScale', 'linear', 'TickLabelInterpreter', 'latex', 'FontSize', 18);
grid on; box on;
title('\textbf{DMD: Deviation from Mean}', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('Lead Time', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Error $-$ Mean', 'Interpreter', 'latex', 'FontSize', 14);

xlim([1, size(er2_all, 2)]);
y_range = max(dev(:)) - min(dev(:));
ylim([min(dev(:)) - 0.1*y_range, max(dev(:)) + 0.1*y_range]);

ax = gca;
ax.YAxis.Exponent = floor(log10(max(abs(dev(:))))); % auto scientific notation

% --- 3. Export Vector Graphics PDF ---
exportgraphics(gcf, 'ensemble_weather_error_lines.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none');


%%   ====================kernel_error_all 

figure
hold on;

% Define color
c_kernel = [0.9290, 0.6940, 0.1250];

x = 0:steps;

% Compute mean, min, max across the ensemble dimension (columns)
m_kernel  = mean(Kx0_er_all, 2)';   % 1 x (steps+1)
lo_kernel = min(Kx0_er_all, [], 2)';
hi_kernel = max(Kx0_er_all, [], 2)';

% Shaded ensemble spread (min-max envelope)
fill([x, fliplr(x)], [lo_kernel, fliplr(hi_kernel)], c_kernel, ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Mean line on top
semilogy(x, m_kernel, 'LineWidth', 3, 'Color', c_kernel, ...
    'DisplayName', 'DualKoop-RKHS')

set(gca, 'YScale', 'log')
grid on; box on;

% Set axis tick mark font size and tick label interpreter
ax = gca;
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';

title('Relative forecast error for kernel function', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Snapshots)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend('interpreter', 'latex', 'fontsize', 16, 'location', 'best')

exportgraphics(gcf, 'ensemble_prediction_kernel_error_all.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none')



%%======================================================
figure('Color', 'w', 'Position', [100, 100, 1400, 450])

c_ked  = [0.8500, 0.3250, 0.0980];
c_kmd  = [0.9290, 0.6940, 0.1250];
c_exct = [0, 0.4470, 0.7410];

x = 0:steps;

% --- Panel 1: DMD  (mean_obs_exact_all, blue, matches original p3) ---
subplot(1,3,1); hold on;
lo_exct = log(min(mean_obs_exact_all, [], 2))';
hi_exct = log(max(mean_obs_exact_all, [], 2))';
fill([x, fliplr(x)], [lo_exct, fliplr(hi_exct)], c_exct, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(x, log(mean(mean_obs_exact_all, 2)), '--', 'linewidth', 2.5, 'color', c_exct);
grid on; box on;
ax = gca; ax.FontSize = 20; ax.TickLabelInterpreter = 'latex';
title('DMD', 'interpreter', 'latex', 'fontsize', 16)
xlabel('Lead Time', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Prediction', 'interpreter', 'latex', 'fontsize', 14)

% --- Panel 2: KeDMD  (mean_obs_kedmd_all, orange, matches original p1) ---
subplot(1,3,2); hold on;
lo_ked = log(min(mean_obs_kedmd_all, [], 2))';
hi_ked = log(max(mean_obs_kedmd_all, [], 2))';
fill([x, fliplr(x)], [lo_ked, fliplr(hi_ked)], c_ked, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(x, log(mean(mean_obs_kedmd_all, 2)), 'linewidth', 2.5, 'color', c_ked);
grid on; box on;
ax = gca; ax.FontSize = 20; ax.TickLabelInterpreter = 'latex';
title('KeDMD', 'interpreter', 'latex', 'fontsize', 16)
xlabel('Lead Time', 'interpreter', 'latex', 'fontsize', 14)

% --- Panel 3: DualKoop-RKHS  (mean_obs_kmd_all, yellow, matches original p2) ---
subplot(1,3,3); hold on;
lo_kmd = log(min(mean_obs_kmd_all, [], 2))';
hi_kmd = log(max(mean_obs_kmd_all, [], 2))';
fill([x, fliplr(x)], [lo_kmd, fliplr(hi_kmd)], c_kmd, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(x, log(mean(mean_obs_kmd_all, 2)), '-.', 'linewidth', 2.5, 'color', c_kmd);
grid on; box on;
ax = gca; ax.FontSize = 20; ax.TickLabelInterpreter = 'latex';
title('DualKoop-RKHS', 'interpreter', 'latex', 'fontsize', 16)
xlabel('Lead Time', 'interpreter', 'latex', 'fontsize', 14)

sgtitle('Forecast Comparison (per-method scale)', 'interpreter', 'latex', 'fontsize', 18)

exportgraphics(gcf, 'ensemble_prediction_mean_panels.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')


%%======================================================
figure('Color', 'w', 'Position', [100, 100, 1400, 450])
c_ked  = [0.8500, 0.3250, 0.0980];
c_kmd  = [0.9290, 0.6940, 0.1250];
c_exct = [0, 0.4470, 0.7410];
x = 0:steps;

% --- Panel 1: DMD ---
subplot(1,3,1); hold on;
lo_exct = log(min(mean_obs_exact_all, [], 2))';
hi_exct = log(max(mean_obs_exact_all, [], 2))';
fill([x, fliplr(x)], [lo_exct, fliplr(hi_exct)], c_exct, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(x, log(mean(mean_obs_exact_all, 2)), '--', 'linewidth', 2.5, 'color', c_exct, 'DisplayName', 'DMD');
grid on; box on;
ax = gca; ax.FontSize = 22; ax.TickLabelInterpreter = 'latex';
%title('DMD', 'interpreter', 'latex', 'fontsize', 22)
xlabel('Lead Time', 'interpreter', 'latex', 'fontsize', 22)
ylabel('Prediction', 'interpreter', 'latex', 'fontsize', 22)
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'best');

% --- Panel 2: KeDMD ---
subplot(1,3,2); hold on;
lo_ked = log(min(mean_obs_kedmd_all, [], 2))';
hi_ked = log(max(mean_obs_kedmd_all, [], 2))';
fill([x, fliplr(x)], [lo_ked, fliplr(hi_ked)], c_ked, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(x, log(mean(mean_obs_kedmd_all, 2)), '-', 'linewidth', 2.5, 'color', c_ked, 'DisplayName', 'KeDMD');
grid on; box on;
ax = gca; ax.FontSize = 22; ax.TickLabelInterpreter = 'latex';
%title('KeDMD', 'interpreter', 'latex', 'fontsize', 22)
xlabel('Lead Time', 'interpreter', 'latex', 'fontsize', 22)
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'best');

% --- Panel 3: DualKoop-RKHS ---
subplot(1,3,3); hold on;
lo_kmd = log(min(mean_obs_kmd_all, [], 2))';
hi_kmd = log(max(mean_obs_kmd_all, [], 2))';
fill([x, fliplr(x)], [lo_kmd, fliplr(hi_kmd)], c_kmd, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(x, log(mean(mean_obs_kmd_all, 2)), '-.', 'linewidth', 2.5, 'color', c_kmd, 'DisplayName', 'DualKoop-RKHS');
grid on; box on;
ax = gca; ax.FontSize = 22; ax.TickLabelInterpreter = 'latex';
%title('DualKoop-RKHS', 'interpreter', 'latex', 'fontsize', 22)
xlabel('Lead Time', 'interpreter', 'latex', 'fontsize', 22)
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'best');

sgtitle('Forecast Comparison (per-method scale)', 'interpreter', 'latex', 'fontsize', 24)
exportgraphics(gcf, 'ensemble_prediction_mean_panels.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% === One palnel ensemble-prediction-mean-shaoow
figure('Color', 'w', 'Position', [100, 100, 900, 600])
hold on;

c_ked  = [0.8500, 0.3250, 0.0980];
c_kmd  = [0.9290, 0.6940, 0.1250];
c_exct = [0, 0.4470, 0.7410];

x = 0:steps;

% --- Shaded bands (min-max across ensemble, log space) ---
% EdgeColor is set (thin) so a very narrow band stays visible as an outline
lo_exct = log(min(mean_obs_exact_all, [], 2))';
hi_exct = log(max(mean_obs_exact_all, [], 2))';
fill([x, fliplr(x)], [lo_exct, fliplr(hi_exct)], c_exct, ...
    'FaceAlpha', 0.25, 'EdgeColor', c_exct, 'LineWidth', 0.5, 'HandleVisibility', 'off');

lo_ked = log(min(mean_obs_kedmd_all, [], 2))';
hi_ked = log(max(mean_obs_kedmd_all, [], 2))';
fill([x, fliplr(x)], [lo_ked, fliplr(hi_ked)], c_ked, ...
    'FaceAlpha', 0.25, 'EdgeColor', c_ked, 'LineWidth', 0.5, 'HandleVisibility', 'off');

lo_kmd = log(min(mean_obs_kmd_all, [], 2))';
hi_kmd = log(max(mean_obs_kmd_all, [], 2))';
fill([x, fliplr(x)], [lo_kmd, fliplr(hi_kmd)], c_kmd, ...
    'FaceAlpha', 0.25, 'EdgeColor', c_kmd, 'LineWidth', 0.5, 'HandleVisibility', 'off');

% --- Mean lines on top ---
p1 = plot(x, log(mean(mean_obs_exact_all, 2)),  '--', 'LineWidth', 2.5, 'Color', c_exct);
p2 = plot(x, log(mean(mean_obs_kedmd_all, 2)),  '-',  'LineWidth', 2.5, 'Color', c_ked);
p3 = plot(x, log(mean(mean_obs_kmd_all, 2)),    '-.', 'LineWidth', 2.5, 'Color', c_kmd);

grid on; box on;
ax = gca; ax.FontSize = 22; ax.TickLabelInterpreter = 'latex';
xlabel('Lead Time', 'Interpreter', 'latex', 'FontSize', 22)
ylabel('Prediction', 'Interpreter', 'latex', 'FontSize', 22)
title('Forecast Comparison', 'Interpreter', 'latex', 'FontSize', 24)
legend([p1 p2 p3], {'DMD', 'KeDMD', 'DualKoop-RKHS'}, ...
    'Interpreter', 'latex', 'FontSize', 18, 'Location', 'best');

exportgraphics(gcf, 'ensemble_prediction_mean_one_panel.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none')

%% Spaltial modes 

%% Specific Location (Kumamoto)
lat_sub = weatherDat2021AUG_ensemble0.dat_lat;  % 140x149
lon_sub = weatherDat2021AUG_ensemble0.dat_lon;  % 140x149
shapefile_path = 'D:\Susuki Lab\Testing_Code\data-weather\Data_250401\ne_10m_coastline\ne_10m_coastline.shp';
S = shaperead(shapefile_path);

targetLat = 32.803; targetLon = 130.708;
dist = sqrt((lat_sub - targetLat).^2 + (lon_sub - targetLon).^2);
[minDist, Kumamoto_Locatidx] = min(dist(:));
[r, c] = ind2sub(size(lat_sub), Kumamoto_Locatidx);
fprintf('Kumamoto is at Index: %d (Row: %d, Col: %d)\n', Kumamoto_Locatidx, r, c);

%% Plot spatial modes
fig = figure('Color','w','Position',[100 100 1200 900]);
colormap(brighten(redblueTecplot(21), -0.55));
mode_indices = [1,3,5,7,9,11,13,15,17];

pp = 1;
for j = mode_indices
    % 1. Eigenvalue, eigenvector, residual
    L   = Lambda_res(j);
    F_j = F_res(:, j);
    r_j = res_verif(j);

    % 2. Reconstruct spatial mode and phase align
    Phi = ((G * F_j) \ (X.')).';
    u   = Phi;
    u   = real(u * exp(1i * mean(angle(u))));

    % 3. Map land vector to 2D grid
    v = nan(size(lon_sub));
    v(nLAND) = u(:);
    mode_map = abs(v);
    mag_L = abs(L);

    % 4. Subplot
    subplot(3, 3, pp);
    p_plot = pcolor(lon_sub, lat_sub, mode_map);
    p_plot.EdgeColor = 'none';
    shading interp;
    axis equal tight;
    set(gca, 'YDir', 'normal');

    c_max = max(mode_map(:), [], 'omitnan');
    c_min = min(mode_map(:), [], 'omitnan');
    if ~isnan(c_min) && ~isnan(c_max) && c_min ~= c_max
        clim([c_min, c_max]);
    end
    colorbar;
    hold on;

    % Coastline
    for kk = 1:length(S)
        plot(S(kk).X, S(kk).Y, 'k-', 'LineWidth', 1.5);
    end

    %% mark Kumamoto
    plot(lon_sub(r,c), lat_sub(r,c), 'p', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
    hold off;

    title(sprintf('$V_{%d}$: $|\\lambda|=%.3f$, Res$=%.3f$', j, mag_L, r_j), ...
        'Interpreter', 'latex', 'FontSize', 10);
    xlabel('Lon ($^\circ$E)', 'Interpreter', 'latex');
    ylabel('Lat ($^\circ$N)', 'Interpreter', 'latex');
    xlim([min(lon_sub(:)), max(lon_sub(:))]);
    ylim([min(lat_sub(:)), max(lat_sub(:))]);

    pp = pp + 1;
end

% sgtitle('Dual Koopman spatial modes $V_j$', ...
%     'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');

exportgraphics(fig, 'top_9_dual_spatial_modes_r2.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none');

%%  Zoom in -shallow- plots 

figure('Color', 'w', 'Position', [100, 100, 900, 650]);
hold on;

% --- 1. Data Preparation ---
N_steps = size(er1_all, 2);
t_vec   = 1:N_steps;
x_poly  = [t_vec, fliplr(t_vec)];

% Compute min-max shaded ranges across ensembles (dim 1)
min1 = min(er1_all, [], 1); max1 = max(er1_all, [], 1);
min2 = min(er2_all, [], 1); max2 = max(er2_all, [], 1);
min3 = min(er3_all, [], 1); max3 = max(er3_all, [], 1);

% Colors
c_dmd  = [0.0000, 0.4470, 0.7410];
c_kedmd  = [0.8500, 0.3250, 0.0980];
c_dual = [0.9290, 0.6940, 0.1250];

% --- 2. Main Plot (All 3 Methods) ---
fill(x_poly, [max2, fliplr(min2)], c_dmd,  'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(x_poly, [max3, fliplr(min3)], c_kedmd,  'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(x_poly, [max1, fliplr(min1)], c_dual, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

plot(mean(er2_all, 1), '--', 'Color', c_dmd,  'LineWidth', 2, 'DisplayName', 'DMD (mean)');
plot(mean(er3_all, 1), '-',  'Color', c_kedmd,  'LineWidth', 2, 'DisplayName', 'KeDMD (mean)');
plot(mean(er1_all, 1), '-.', 'Color', c_dual, 'LineWidth', 2, 'DisplayName', 'DualKoop-RKHS (mean)');

% Main Axes Formatting
set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 25);
grid on; box on;
title('Relative Forecast Errors Comparison', 'Interpreter', 'latex', 'FontSize', 18);
xlabel('Lead Time (Snapshots)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Relative Forecast Error', 'Interpreter', 'latex', 'FontSize', 18);
legend('Interpreter', 'latex', 'FontSize', 15, 'Location', 'northwest');

% --- 3. Inset Zoom Plot (KeDMD vs DualKoop-RKHS) ---
% Position: [left, bottom, width, height]
axes('Position', [0.50, 0.22, 0.38, 0.38]);
hold on;

% Plot shaded areas for KeDMD and DualKoop
fill(x_poly, [max3, fliplr(min3)], c_kedmd,  'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(x_poly, [max1, fliplr(min1)], c_dual, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Plot mean lines for KeDMD and DualKoop
plot(mean(er3_all, 1), '-',  'Color', c_kedmd,  'LineWidth', 2);
plot(mean(er1_all, 1), '-.', 'Color', c_dual, 'LineWidth', 2);

% Inset Formatting
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
grid on; box on;
title('\textbf{Zoomed: KeDMD vs DualKoop}', 'Interpreter', 'latex', 'FontSize', 14);

% Automatically adjust inset Y-limits based on KeDMD and DualKoop values
y_min_zoom = min([min1, min3]);
y_max_zoom = max([max1, max3]);
xlim([1, N_steps]);
ylim([y_min_zoom * 0.9, y_max_zoom * 1.1]);

% --- 4. Export Vector Graphics PDF ---
exportgraphics(gcf, 'ensemble_shallow_weather_error_1.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none');


%%============================== Lines plots 

figure('Color', 'w', 'Position', [100, 100, 900, 650]);
hold on;

% Define colors
c_dmd  = [0.4, 0.4, 0.4];
c_kedmd  = [0.8500, 0.3250, 0.0980];
c_dual = [0.2, 0.4, 0.8];
% c_dmd  = [0.0000, 0.4470, 0.7410];
% c_kedmd  = [0.8500, 0.3250, 0.0980];
% c_dual = [0.9290, 0.6940, 0.1250];

% --- 1. Main Plot (All Ensembles + Means) ---
M = size(er1_all, 1);
for m = 1:M
    plot(er2_all(m,:), 'Color', [c_dmd, 0.15],  'LineWidth', 1, 'HandleVisibility', 'off')
    plot(er3_all(m,:), 'Color', [c_kedmd, 0.15],  'LineWidth', 1, 'HandleVisibility', 'off')
    plot(er1_all(m,:), 'Color', [c_dual, 0.15], 'LineWidth', 1, 'HandleVisibility', 'off')
end

plot(mean(er2_all, 1), '--', 'Color', c_dmd,  'LineWidth', 2, 'DisplayName', 'DMD (mean)')
plot(mean(er3_all, 1), '-',  'Color', c_kedmd,  'LineWidth', 2, 'DisplayName', 'KeDMD (mean)')
plot(mean(er1_all, 1), '-.', 'Color', c_dual, 'LineWidth', 2, 'DisplayName', 'DualKoop-RKHS (mean)')

% Main Axes Formatting
set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 25);
grid on; box on;
title('Relative Forecast Errors Comparison', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('Lead Time (Snapshots)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Relative Forecast Error', 'Interpreter', 'latex', 'FontSize', 20);
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'northwest');


% --- 2. Inset Zoom Plot (KeDMD vs DualKoop-RKHS in Log Scale) ---
% Position: [left, bottom, width, height]
axes('Position', [0.48, 0.22, 0.40, 0.38]);
hold on;

% Plot ensemble lines for KeDMD and DualKoop only
for m = 1:M
    plot(er3_all(m,:), 'Color', [c_kedmd, 0.2],  'LineWidth', 0.8, 'HandleVisibility', 'off');
    plot(er1_all(m,:), 'Color', [c_dual, 0.2], 'LineWidth', 0.8, 'HandleVisibility', 'off');
end

% Plot mean lines
plot(mean(er3_all, 1), '-',  'Color', c_kedmd,  'LineWidth', 2);
plot(mean(er1_all, 1), '-.', 'Color', c_dual, 'LineWidth', 2);

% Inset Formatting (Log Scale)
set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 18);
grid on; box on;
title('\textbf{Zoomed: KeDMD vs DualKoop}', 'Interpreter', 'latex', 'FontSize', 18);

% Auto-tight y-limits around KeDMD and DualKoop data
y_min_zoom = min([er1_all(:); er3_all(:)]);
y_max_zoom = max([er1_all(:); er3_all(:)]);
xlim([1, size(er1_all, 2)]);
ylim([y_min_zoom * 0.9, y_max_zoom * 1.1]);


% --- 3. Export Vector Graphics PDF ---
exportgraphics(gcf, 'ensemble_weather_error_lines.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none');


%%%% bar 

%% Bar plot: Relatev error


%% ==============Bar ZOOM-------------
c_dmd  = [0.4, 0.4, 0.4];
c_kedmd  = [0.8500, 0.3250, 0.0980];
c_dual = [0.2, 0.4, 0.8];
c_spec= [0.2, 0.4, 0.8];
figure('Color', 'w', 'Position', [100, 100, 800, 600]);

% --- 1. Main Bar Plot ---
b = bar(1:3, rmse_values, 'FaceColor', 'flat', 'BarWidth', 0.5);
b.CData(1,:) = c_dmd;
b.CData(2,:) = c_kedmd;
b.CData(3,:) = c_spec;

% Apply LaTeX interpreter to X-tick labels and axes
set(gca, 'TickLabelInterpreter', 'latex', ...
    'XTick', 1:3, ...
    'XTickLabel', {'\textrm{DMD}', '\textrm{KeDMD}', '\textrm{DualKoop}'});
grid on; box on;
ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 18);
title('Forecast Error Comparison', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 22);


% --- 2. Inset Plot for Small Values (Shifted Slightly Down) ---
% Position: [left, bottom, width, height]
axes('Position', [0.50, 0.32, 0.38, 0.38]); % Changed bottom from 0.45 to 0.32
b_inset = bar(2:3, rmse_values(2:3), 'FaceColor', 'flat', 'BarWidth', 0.5);
b_inset.CData(1,:) = c_kedmd;
b_inset.CData(2,:) = c_spec;

% Apply LaTeX interpreter to inset plot
set(gca, 'TickLabelInterpreter', 'latex', ...
    'XTick', 2:3, ...
    'XTickLabel', {'\textrm{KeDMD}', '\textrm{DualKoop}'});
title('\textbf{Zoomed (KeDMD vs DualKoop)}', 'Interpreter', 'latex', 'FontSize', 11);
grid on; box on;
set(gca, 'FontSize', 14);


% --- 3. Export Vector Graphics PDF ---
exportgraphics(gcf, 'ensemble_overall_rmse_bar_1.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'none');


%% Bar plot: Relatev error
% % --- 1. Compute overall RMSE across ALL time steps and ALL ensembles ---
% % (er(:).^2 squares all elements, mean computes average, sqrt gives RMSE)
% rmse_DMD      = sqrt(mean(er2_all(:).^2));
% rmse_kEDMD    = sqrt(mean(er3_all(:).^2));
% rmse_SpecRKHS = sqrt(mean(er1_all(:).^2));
% 
% rmse_values = [rmse_DMD; rmse_kEDMD; rmse_SpecRKHS];
% 
% % --- 2. Color definitions (Matching your line plot) ---
% c_dmd   = [0.0000 0.4470 0.7410]; % Blue
% c_kedmd = [0.8500 0.3250 0.0980]; % Terracotta / Red
% c_spec  = [0.9290 0.6940 0.1250]; % Orange
% 
% % --- 3. Plot Bar Chart ---
% figure('Color', 'w');
% b = bar(1:3, rmse_values, 'FaceColor', 'flat', 'BarWidth', 0.5);
% b.CData(1,:) = c_dmd;
% b.CData(2,:) = c_kedmd;
% b.CData(3,:) = c_spec;
% 
% % --- 4. Log Scale & Formatting ---
% set(gca, 'YScale', 'log');
% set(gca, 'XTick', 1:3, 'XTickLabel', {'DMD', 'KeDMD', 'DualKoop'});
% grid on; box on;
% ylabel('Overall RMSE', 'interpreter', 'latex', 'fontsize', 18);
% title('Overall Forecast Error Comparison (RMSE)', 'interpreter', 'latex', 'fontsize', 18);
% set(gca, 'FontSize', 14);
% 
% % Export vector graphics PDF
% exportgraphics(gcf, 'ensemble_overall_rmse_bar.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
