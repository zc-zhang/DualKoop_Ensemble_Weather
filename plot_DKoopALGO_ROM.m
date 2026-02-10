% Low-dimensional system modeling.

%% ========================================================================
%% STEP 1: Build Dual Koopman Operator (You already have this)
%% ========================================================================
[G, K_star, L, lambda, W, C_ref, PSI_ref, ~, ~, ~, UU_ref] = dualKoop_Algom( ...
    Xa_mean, Ya_mean, ...
    'type', 'Linear', ...
    'N', 300, ...
    'cut_off', 1e-10 ...
);

%% ========================================================================
%% STEP 2: SELECT LOW-DIMENSIONAL MODES (Model Reduction!)
%% ========================================================================
% Choose number of modes to keep
r = 70;  % Low-dimensional approximation (r << N)

% Extract dominant eigenvalues and modes
[~, idx_sorted] = sort(abs(lambda), 'descend');
lambda_r = lambda(idx_sorted(1:r));          % Top r eigenvalues
C_r = C_ref(:, idx_sorted(1:r));             % Top r dual Koopman modes (N × r)

% Dynamics matrix (diagonal)
Lambda_r = diag(lambda_r);  % r × r

fprintf('\n=== Low-Dimensional Model ===\n');
fprintf('Original dimension: %d\n', size(Xa_mean, 1));
fprintf('Reduced dimension: r = %d\n', r);
fprintf('Compression ratio: %.2f%%\n', 100*r/size(Xa_mean, 1));

%% ========================================================================
%% STEP 3: COMPUTE REDUCED COORDINATES FOR ALL ENSEMBLES
%% ========================================================================
% For each ensemble, compute α^(i)(k) ∈ ℝ^r
alpha_ens = cell(M, 1);    % Reduced coordinates
Y_recon_ens = cell(M, 1);  % Reconstructed outputs

fprintf('\nComputing reduced coordinates for %d ensembles...\n', M);

parfor m = 1:M
    % Evaluate eigenfunctions on this ensemble
    [~, ~, ~, ~, ~, ~, PSI_m] = dualKoop_Algom( ...
        Xa_mean, Ya_mean, ...
        'Xb', Xa_ens{m}, ...
        'UU', UU_ref ...
    );
    
    % Extract only the r dominant eigenfunctions
    PSI_m_r = PSI_m(:, idx_sorted(1:r));  % T × r
    
    % These ARE the reduced coordinates!
    % α^(i)(k) = [φ_1(x_k^(i)), φ_2(x_k^(i)), ..., φ_r(x_k^(i))]^T
    alpha_ens{m} = PSI_m_r;  % T × r matrix (each row is α at time k)
    
    % Reconstruct output (if you have observable data)
    % This would require the B matrix: y_k ≈ B * α_k
    % For now, we'll focus on the reduced state dynamics
end

%% ========================================================================
%% STEP 4: COMPUTE ENSEMBLE MEAN IN REDUCED SPACE
%% ========================================================================
% Align time dimensions
min_T = min(cellfun(@(x) size(x, 1), alpha_ens));
alpha_aligned = cellfun(@(a) a(1:min_T, :), alpha_ens, 'UniformOutput', false);

% Stack into 3D array: T × r × M
alpha_all = cat(3, alpha_aligned{:});  % T × r × M

% Ensemble mean reduced coordinates
alpha_mean = mean(alpha_all, 3);  % T × r

% Ensemble standard deviation
alpha_std = std(alpha_all, 0, 3);  % T × r

fprintf('\n=== Ensemble Statistics in Reduced Space ===\n');
fprintf('Time steps: %d\n', min_T);
fprintf('Reduced dimension: %d\n', r);
fprintf('Number of ensembles: %d\n', M);

%% ========================================================================
%% STEP 5: DEVIATION ANALYSIS IN REDUCED SPACE
%% ========================================================================
% Compute deviation from ensemble mean
deviation_norm = zeros(M, min_T);

for m = 1:M
    for k = 1:min_T
        alpha_dev = alpha_all(k, :, m) - alpha_mean(k, :);
        deviation_norm(m, k) = norm(alpha_dev);
    end
end

%% ========================================================================
%% STEP 6: VISUALIZATION - LOW-DIMENSIONAL DYNAMICS
%% ========================================================================

%% Plot 1: Reduced Coordinates Over Time
figure('Position', [100, 100, 1400, 800]);

% Plot first 3 modes for selected ensembles
ensemble_subset = [1, 25, 50, 75, 100];
colors = lines(length(ensemble_subset));

for j = 1:min(3, r)  % Plot first 3 modes
    subplot(2, 3, j);
    hold on;
    
    for i = 1:length(ensemble_subset)
        m = ensemble_subset(i);
        plot(real(alpha_all(:, j, m)), 'Color', colors(i, :), ...
            'LineWidth', 1.5, 'DisplayName', sprintf('Ens %03d', m-1));
    end
    
    % Plot ensemble mean
    plot(real(alpha_mean(:, j)), 'k--', 'LineWidth', 2.5, ...
        'DisplayName', 'Mean');
    
    xlabel('Time index k');
    ylabel(sprintf('\\alpha_%d^{(i)}(k)', j));
    title(sprintf('Mode %d (\\lambda_%d = %.4f)', j, j, lambda_r(j)));
    legend('Location', 'best');
    grid on;
    hold off;
end

%% Plot 2: Ensemble Spread Over Time
subplot(2, 3, 4);
hold on;
for j = 1:min(3, r)
    plot(alpha_std(:, j), 'LineWidth', 2, ...
        'DisplayName', sprintf('Mode %d', j));
end
xlabel('Time index k');
ylabel('Ensemble std dev');
title('Ensemble Uncertainty in Reduced Space');
legend('Location', 'best');
grid on;

%% Plot 3: Total Deviation from Mean
subplot(2, 3, 5);
imagesc(deviation_norm');
colorbar;
xlabel('Ensemble member m');
ylabel('Time index k');
title('||α^{(i)}(k) - \barα(k)||');
set(gca, 'YDir', 'normal');

%% Plot 4: Time-averaged deviation per ensemble
subplot(2, 3, 6);
mean_deviation = mean(deviation_norm, 2);
bar(0:M-1, mean_deviation);
xlabel('Ensemble member m');
ylabel('Time-averaged deviation');
title('Ensemble Member Deviation from Mean');
grid on;

%% ========================================================================
%% STEP 7: MODE CONTRIBUTION ANALYSIS
%% ========================================================================
figure('Position', [100, 100, 1200, 500]);

% Subplot 1: Eigenvalue spectrum
subplot(1, 3, 1);
semilogy(1:r, abs(lambda_r), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Mode index j');
ylabel('|λ_j|');
title('Retained Eigenvalue Spectrum');
grid on;

% Subplot 2: Variance explained by each mode
mode_variance = zeros(r, 1);
for j = 1:r
    mode_variance(j) = var(alpha_all(:, j, :), 0, 'all');
end
mode_variance = mode_variance / sum(mode_variance);  % Normalize

subplot(1, 3, 2);
bar(1:r, 100*mode_variance);
xlabel('Mode index j');
ylabel('% Variance explained');
title('Mode Contribution to Total Variance');
grid on;

% Subplot 3: Cumulative variance
subplot(1, 3, 3);
plot(1:r, 100*cumsum(mode_variance), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of modes');
ylabel('% Cumulative variance');
title('Cumulative Variance Explained');
grid on;
ylim([0, 105]);

%% ========================================================================
%% STEP 8: DEVIATION ERROR BOUNDS
%% ========================================================================
figure('Position', [100, 100, 1400, 600]);

% Compute reconstruction error for different r values
%r_values = [5, 10, 15, 20, 30, 40, 50, 60, 70];
r_values = [5, 10, 15, 20, 30, 50];
reconstruction_error = zeros(length(r_values), M);

for i_r = 1:length(r_values)
    r_test = r_values(i_r);
    if r_test > length(lambda)
        continue;
    end
    
    % Use only first r_test modes
    for m = 1:M
        alpha_full = alpha_aligned{m};  % Full reduced coords
        alpha_truncated = alpha_full(:, 1:r_test);
        
        % Error is norm of discarded modes
        if size(alpha_full, 2) > r_test
            alpha_residual = alpha_full(:, r_test+1:end);
            reconstruction_error(i_r, m) = norm(alpha_residual, 'fro') / norm(alpha_full, 'fro');
        end
    end
end

% Plot error vs number of modes
subplot(1, 2, 1);
boxplot(reconstruction_error', r_values);
xlabel('Number of modes r');
ylabel('Relative reconstruction error');
title('Truncation Error vs Model Dimension');
grid on;
set(gca, 'YScale', 'log');

% Plot error bound envelope
subplot(1, 2, 2);
mean_error = mean(reconstruction_error, 2);
std_error = std(reconstruction_error, 0, 2);

errorbar(r_values, mean_error, std_error, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of modes r');
ylabel('Mean ± std reconstruction error');
title('Error Bounds for Model Reduction');
grid on;
set(gca, 'YScale', 'log');

%% ========================================================================
%% STEP 9: STATISTICS OUTPUT
%% ========================================================================
fprintf('\n=== Low-Dimensional Model Performance ===\n');
fprintf('Reduced dimension r = %d\n', r);
fprintf('\nMode stability:\n');
for j = 1:min(5, r)
    stability = abs(lambda_r(j));
    if stability < 1
        trend = 'decaying';
    elseif stability > 1
        trend = 'growing';
    else
        trend = 'neutral';
    end
    fprintf('  Mode %d: |λ_%d| = %.4f (%s)\n', j, j, stability, trend);
end

fprintf('\nEnsemble spread (time-averaged):\n');
for j = 1:min(5, r)
    fprintf('  Mode %d: std = %.4e\n', j, mean(alpha_std(:, j)));
end

fprintf('\nDeviation from ensemble mean:\n');
fprintf('  Mean deviation: %.4e\n', mean(mean_deviation));
fprintf('  Max deviation:  %.4e (ensemble %d)\n', ...
    max(mean_deviation), find(mean_deviation == max(mean_deviation)) - 1);
fprintf('  Min deviation:  %.4e (ensemble %d)\n', ...
    min(mean_deviation), find(mean_deviation == min(mean_deviation)) - 1);

%% ========================================================================
%% STEP 10: IDENTIFY "BEST" AND "WORST" ENSEMBLE MEMBERS
%% ========================================================================
[~, idx_best] = min(mean_deviation);
[~, idx_worst] = max(mean_deviation);

fprintf('\n=== Ensemble Member Classification ===\n');
fprintf('Most representative (closest to mean): enPREC%03d\n', idx_best - 1);
fprintf('Most extreme (farthest from mean):     enPREC%03d\n', idx_worst - 1);

% Plot comparison
figure('Position', [100, 100, 1200, 400]);
for j = 1:min(3, r)
    subplot(1, 3, j);
    hold on;
    plot(real(alpha_all(:, j, idx_best)), 'b-', 'LineWidth', 2, ...
        'DisplayName', sprintf('Best (Ens %03d)', idx_best-1));
    plot(real(alpha_all(:, j, idx_worst)), 'r-', 'LineWidth', 2, ...
        'DisplayName', sprintf('Worst (Ens %03d)', idx_worst-1));
    plot(real(alpha_mean(:, j)), 'k--', 'LineWidth', 2, ...
        'DisplayName', 'Mean');
    xlabel('Time k');
    ylabel(sprintf('\\alpha_%d(k)', j));
    title(sprintf('Mode %d Comparison', j));
    legend('Location', 'best');
    grid on;
end

fprintf('\nAnalysis complete!\n');