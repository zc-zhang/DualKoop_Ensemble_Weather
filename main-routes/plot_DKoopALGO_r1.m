% plot_DKoopALGO_r1.m
% plot_DKoopALGO_corrected.m
% Dual Koopman Eigenfunction Ensemble Analysis

% Mathematical Framework:
% ε_t^(i) = x_t^(i) - x̄_t = Σ_j λ_j^t [φ_j(x_0^(i)) - φ_j(x̄_0)] c_j
%
% where:
%   - x_t^(i): state of ensemble i at time t
%   - x̄_t: mean ensemble state
%   - φ_j: dual Koopman eigenfunctions
%   - λ_j: Koopman eigenvalues
%   - c_j: dual Koopman modes


% Dual Koopman Deviation Dynamics Analysis
% 
% Mathematical Framework:
% ε_{t,ℓ}^(i) = y_{t,ℓ}^(i) - ȳ_{t,ℓ}
%             = Σ_j λ_j^t [φ_j(x_0^(i)) - φ_j(x̄_0)] ⟨ζ_j, h_ℓ⟩
%
% Data Structure:
%   - Xraw: 20860 × 73 (time × space)
%   - Xa_ens{m}: 20860 × 71 (remove boundary spatial points)
%   - enPREC100 = mean trajectory
% =========================================================================

%clear; clc;

%% ============================================================
%% STEP 1: Load Ensemble Data
%% ============================================================
%fprintf('=== DUAL KOOPMAN DEVIATION ANALYSIS ===\n\n');

data_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Ensemble_PREC';
M = 100;  % Ensemble members: enPREC000 to enPREC099
Xa_ens = cell(M, 1);  
Ya_ens = cell(M, 1);

fprintf('Loading %d ensemble members (enPREC000-099)...\n', M);
for m = 0:M-1
    fname = sprintf('enPREC%03d.mat', m);
    S = load(fullfile(data_dir, fname));
    Xraw = S.vectorized_PREC;    % 20860 × 73 (time × space)
    
    % Remove boundary spatial points (columns 1 and 73)
    % Keep columns 2 to 72 (71 interior spatial points)
    Xa_ens{m+1} = Xraw(:, 2:end-1);  % 20860 × 71: x^(i)(t)
    Ya_ens{m+1} = Xraw(:, 3:end);    % 20860 × 72: shifted spatial pattern
end

%% Load MEAN trajectory (enPREC100)
fprintf('Loading mean trajectory (enPREC100)...\n');
S_mean = load(fullfile(data_dir, 'enPREC100.mat'));
Xraw_mean = S_mean.vectorized_PREC;  % 20860 × 73
Xa_mean = Xraw_mean(:, 2:end-1);     % 20860 × 71: x̄(t)
Ya_mean = Xraw_mean(:, 3:end);       % 20860 × 72

fprintf('Data loaded successfully.\n');
fprintf('  Ensemble members: %d\n', M);
fprintf('  Mean trajectory: enPREC100\n');
[nT, nSpace] = size(Xa_mean);
fprintf('  Data format: %d time steps × %d spatial points\n\n', nT, nSpace);

%% Verify dimensions
fprintf('Dimension verification:\n');
fprintf('  Xa_mean: [%d, %d]\n', size(Xa_mean));
fprintf('  Ya_mean: [%d, %d]\n', size(Ya_mean));
fprintf('  Xa_ens{1}: [%d, %d]\n', size(Xa_ens{1}));
fprintf('  Ya_ens{1}: [%d, %d]\n\n', size(Ya_ens{1}));

%% ============================================================
%% STEP 2: Build Dual Koopman Operator on Mean Trajectory
%% Kernel dictionary: {k(·, x̄_i)} based on mean dynamics
%% ============================================================
fprintf('Building dual Koopman operator on mean trajectory (enPREC100)...\n');

[G, A, K_star, lambda, W, C_mean, PSI_mean, PSI_y_mean, ~, kernel_f, UU_mean] = ...
    dualKoop_Algom(Xa_mean, Ya_mean, ...
    'type', 'Gaussian', ...
    'N', 300, ...
    'cut_off', 1e-10, ...
    'Xb', Xa_mean);  % Evaluate on mean trajectory

fprintf('Operator built successfully!\n');
fprintf('  G (Gram matrix): [%d, %d]\n', size(G));
fprintf('  K_star (dual Koopman): [%d, %d]\n', size(K_star));
fprintf('  Number of eigenvalues: %d\n', length(lambda));
fprintf('  C_mean (dual modes): [%d, %d]\n', size(C_mean));
fprintf('  PSI_mean (eigenfunctions): [%d, %d]\n', size(PSI_mean));

% Verify
if isempty(lambda)
    error('Lambda is empty! Check output order: [G, A, K_star, lambda, ...]');
end
if isempty(PSI_mean)
    error('PSI_mean is empty! Ensure Xb=Xa_mean is passed.');
end

fprintf('\n  Top 10 eigenvalues:\n');
for j = 1:min(10, length(lambda))
    fprintf('    λ_%d = %.6f + %.6fi (|λ| = %.6f)\n', ...
        j, real(lambda(j)), imag(lambda(j)), abs(lambda(j)));
end
fprintf('\n');

%% ============================================================
%% STEP 3: Evaluate Eigenfunctions on Each Ensemble
%% Compute: φ_j(x^(i)(t)) for i = 1,...,100
%% ============================================================
fprintf('Evaluating eigenfunctions on %d ensemble trajectories...\n', M);
PSI_ens = cell(M, 1);

parfor m = 1:M
    % Use SAME kernel dictionary (Xa_mean, Ya_mean)
    % Use SAME eigenvectors (UU_mean)
    % Evaluate on ensemble m: Xa_ens{m}
    [~, ~, ~, ~, ~, ~, PSI_ens{m}] = dualKoop_Algom( ...
        Xa_mean, Ya_mean, ...  % Reference (mean) for kernel
        'Xb', Xa_ens{m}, ...   % Evaluate on ensemble i
        'UU', UU_mean ...      % Reuse eigenvectors
    );
    
    if mod(m, 20) == 0
        fprintf('  Processed %d/%d ensembles\n', m, M);
    end
end

fprintf('Eigenfunction evaluation complete.\n');

% Verify
if any(cellfun(@isempty, PSI_ens))
    idx_empty = find(cellfun(@isempty, PSI_ens), 1);
    error('PSI_ens{%d} is empty! Check evaluation step.', idx_empty);
end
fprintf('  PSI_ens{1}: [%d, %d]\n\n', size(PSI_ens{1}));

%% ============================================================
%% STEP 4: Handle Size Differences (if any)
%% ============================================================
sizes = cellfun(@(x) size(x, 1), PSI_ens);
min_size = min([sizes; size(PSI_mean, 1)]);

fprintf('Aligning to minimum size: %d\n', min_size);
PSI_aligned = cellfun(@(P) P(1:min_size, :), PSI_ens, 'UniformOutput', false);
PSI_mean_aligned = PSI_mean(1:min_size, :);

%% ============================================================
%% STEP 5: Eigenfunction Deviation Analysis
%% Δφ_j^(i)(t) = φ_j(x^(i)(t)) - φ_j(x̄(t))
%% ============================================================
fprintf('Computing eigenfunction deviations...\n');

% Select modes for analysis (skip mode 1, often trivial)
eig_indices = [2, 4, 6, 8];
K = length(eig_indices);

% Build eigenfunction matrices
% PHI_matrix(k, t, m) = φ_j(x^(m)(t))
PHI_matrix = zeros(K, min_size, M);
PHI_mean_matrix = zeros(K, min_size);

for k = 1:K
    eig_id = eig_indices(k);
    
    % Mean eigenfunction trajectory: φ_j(x̄(t))
    phi_mean_k = PSI_mean_aligned(:, eig_id);
    
    % Handle complex eigenfunctions
    if ~isreal(phi_mean_k)
        fprintf('  Warning: Mode φ_%d is complex. Taking real part.\n', eig_id);
        phi_mean_k = real(phi_mean_k);
    end
    
    PHI_mean_matrix(k, :) = phi_mean_k';
    
    % Each ensemble: φ_j(x^(i)(t))
    for m = 1:M
        phi_m_k = PSI_aligned{m}(:, eig_id);
        if ~isreal(phi_m_k)
            phi_m_k = real(phi_m_k);
        end
        PHI_matrix(k, :, m) = phi_m_k';
    end
end

% Normalize eigenfunctions
for k = 1:K
    % Normalize mean
    norm_mean = norm(PHI_mean_matrix(k, :));
    if norm_mean > eps
        PHI_mean_matrix(k, :) = PHI_mean_matrix(k, :) / norm_mean;
    end
    
    % Normalize each ensemble
    for m = 1:M
        norm_m = norm(squeeze(PHI_matrix(k, :, m)));
        if norm_m > eps
            PHI_matrix(k, :, m) = PHI_matrix(k, :, m) / norm_m;
        end
    end
end

% Align signs to mean
for k = 1:K
    phi_mean_k = squeeze(PHI_mean_matrix(k, :));
    for m = 1:M
        phi_m_k = squeeze(PHI_matrix(k, :, m));
        if dot(phi_m_k, phi_mean_k) < 0
            PHI_matrix(k, :, m) = -PHI_matrix(k, :, m);
        end
    end
end

% Recompute mean after alignment
for k = 1:K
    PHI_mean_matrix(k, :) = mean(PHI_matrix(k, :, :), 3);
    norm_mean = norm(PHI_mean_matrix(k, :));
    if norm_mean > eps
        PHI_mean_matrix(k, :) = PHI_mean_matrix(k, :) / norm_mean;
    end
end

% Deviation: Δφ_j^(i)(t) = φ_j(x^(i)(t)) - φ_j(x̄(t))
DELTA_PHI = zeros(K, min_size, M);
for m = 1:M
    for k = 1:K
        DELTA_PHI(k, :, m) = PHI_matrix(k, :, m) - PHI_mean_matrix(k, :);
    end
end

% Temporal norm: ||Δφ_j^(i)||_2
E_phi = zeros(K, M);
for k = 1:K
    for m = 1:M
        E_phi(k, m) = norm(squeeze(DELTA_PHI(k, :, m)));
    end
end

% Normalized deviation: Ẽ_j^(i) = ||Δφ_j^(i)||_2 / ||φ_j(x̄)||_2
E_phi_normalized = zeros(K, M);
for k = 1:K
    phi_mean_norm = norm(PHI_mean_matrix(k, :));
    if phi_mean_norm > eps
        E_phi_normalized(k, :) = E_phi(k, :) / phi_mean_norm;
    end
end

% Total deviation energy: E_total^(i) = √(Σ_j (Ẽ_j^(i))²)
E_total = sqrt(sum(E_phi_normalized.^2, 1));

% Mode contribution: γ_j^(i) = (Ẽ_j^(i))² / Σ_k (Ẽ_k^(i))²
gamma = zeros(K, M);
for m = 1:M
    sum_sq = sum(E_phi_normalized(:, m).^2);
    if sum_sq > eps
        gamma(:, m) = E_phi_normalized(:, m).^2 / sum_sq * 100;
    end
end

fprintf('Eigenfunction deviation analysis complete.\n\n');

%% ============================================================
%% STEP 6: State-Space Deviation
%% ε^(i)(t) = x^(i)(t) - x̄(t)
%% ============================================================
fprintf('Computing state-space deviations...\n');

EPSILON_state = zeros(min_size, nSpace, M);
for m = 1:M
    % Align sizes
    Xa_m = Xa_ens{m}(1:min_size, :);
    Xa_mean_cut = Xa_mean(1:min_size, :);
    
    EPSILON_state(:, :, m) = Xa_m - Xa_mean_cut;
end

% Frobenius norm
epsilon_norm = zeros(M, 1);
for m = 1:M
    epsilon_norm(m) = norm(squeeze(EPSILON_state(:, :, m)), 'fro');
end

fprintf('State-space deviations computed.\n\n');

%% ============================================================
%% STEP 7: Pairwise Deviation Matrix
%% D_k^(i,j) = ||φ_k(x^(i)) - φ_k(x^(j))|| / ||φ_k(x^(j))||
%% ============================================================
fprintf('Computing pairwise deviations...\n');

D_pairwise = zeros(K, M, M);
for k = 1:K
    for i = 1:M
        for j = 1:M
            if i == j
                D_pairwise(k, i, j) = 0;
            else
                phi_i = squeeze(PHI_matrix(k, :, i));
                phi_j = squeeze(PHI_matrix(k, :, j));
                D_pairwise(k, i, j) = norm(phi_i - phi_j) / norm(phi_j);
            end
        end
    end
end

% Average pairwise distance for each mode
D_pairwise_mean = zeros(K, 1);
for k = 1:K
    D_k = squeeze(D_pairwise(k, :, :));
    upper_tri = triu(D_k, 1);
    vals = upper_tri(upper_tri > 0);
    D_pairwise_mean(k) = mean(vals);
end

fprintf('Pairwise deviations computed.\n\n');

%% ============================================================
%% STEP 8: Visualization
% plot dual Koopman eigenvalues -mean case
figure('Position', [100, 100, 800, 600]);
scatter(real(lambda), imag(lambda), 80, 'bo', 'filled', 'MarkerEdgeColor', 'b');
hold on;
% Unit circle for reference
theta = linspace(0, 2*pi, 100);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 2);
hold off;

xlabel('Re(\lambda)', 'FontSize', 12);
ylabel('Im(\lambda)', 'FontSize', 12);
title('Dual Koopman Eigenvalue (Mean Ensemble)', 'FontSize', 14);
grid on;
axis equal;
box
%legend('Eigenvalues', 'Unit Circle', 'Location', 'best');
%% ============================================================
fprintf('Generating visualizations...\n');

figure('Position', [50, 50, 1800, 1000]);

%% Plot 1: Eigenvalue Spectrum
subplot(3, 4, 1);
stem(1:min(71, length(lambda)), abs(lambda(1:min(71, length(lambda)))), ...
    'o-', 'LineWidth', 1.5, 'Color', [0.2, 0.4, 0.7]);
hold on;
for k = 1:K
    plot(eig_indices(k), abs(lambda(eig_indices(k))), 'ro', ...
        'MarkerSize', 10, 'LineWidth', 2);
end
xlabel('Mode index j');
ylabel('|\lambda_j|');
title('Magnitude of Dual Koopman Spectra');
grid on;
%legend('All modes', 'Selected', 'Location', 'best');

%% Plot 2-5: Eigenfunction Deviations Δφ_j^(i)(t)
for k = 1:4
    subplot(3, 4, 1+k);
    eig_id = eig_indices(k);
    
    hold on;
    sample_ens = [1, 25, 50, 75, 100];
    colors = lines(length(sample_ens));
    
    for idx = 1:length(sample_ens)
        m = sample_ens(idx);
        delta_phi = squeeze(DELTA_PHI(k, :, m));
        plot(delta_phi, 'Color', colors(idx, :), 'LineWidth', 1.3, ...
            'DisplayName', sprintf('enPREC%03d', m-1));
    end
    yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    xlabel('Time t');
    ylabel(['\Delta\phi_{', num2str(eig_id), '}^{(i)}(t)']);
    title(sprintf('Mode %d: \\lambda = %.4f', eig_id, lambda(eig_id)));
    legend('Location', 'best', 'FontSize', 8);
    grid on;
end

%% Plot 6: Normalized Deviation Heatmap
subplot(3, 4, 6);
imagesc(1:M, eig_indices, E_phi_normalized);
colorbar;
xlabel('Ensemble i (enPREC000-099)');
ylabel('Mode j');
title('Normalized Deviation \tilde{E}_j^{(i)}');
set(gca, 'YDir', 'normal');
colormap('jet');

%% Plot 7: Total Deviation Energy
subplot(3, 4, 7);
bar(0:M-1, E_total, 'FaceColor', [0.3, 0.6, 0.8]);
xlabel('Ensemble i');
ylabel('E_{total}^{(i)}');
title('Total Deviation Energy');
grid on;
xlim([-1, M]);

%% Plot 8: Mode Contribution (stacked)
subplot(3, 4, 8);
bar(0:M-1, gamma', 'stacked');
xlabel('Ensemble i');
ylabel('Contribution (%)');
title('Mode Contribution \gamma_j^{(i)}');
legend(arrayfun(@(x) sprintf('\\phi_%d', x), eig_indices, 'UniformOutput', false), ...
    'Location', 'eastoutside', 'FontSize', 8);
ylim([0, 100]);
grid on;

%% Plot 9: Correlation - State vs Eigenfunction Deviation
subplot(3, 4, 9);
scatter(epsilon_norm, E_total, 60, 'filled', 'MarkerFaceAlpha', 0.7);
xlabel('State deviation ||\epsilon^{(i)}||_F');
ylabel('Eigenfunction deviation E_{total}^{(i)}');
title('Correlation: State vs Eigenfunction');
grid on;
[rho, pval] = corr(epsilon_norm, E_total');
text(0.05, 0.95, sprintf('\\rho = %.3f\np = %.2e', rho, pval), ...
    'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white', ...
    'EdgeColor', 'k');

%% Plot 10: Pairwise Discrepancy Matrix (Mode 2)
subplot(3, 4, 10);
imagesc(0:M-1, 0:M-1, squeeze(D_pairwise(1, :, :)));
colorbar;
xlabel('Ensemble j');
ylabel('Ensemble i');
title(sprintf('Pairwise D_{%d}^{(i,j)}', eig_indices(1)));
axis square;
set(gca, 'YDir', 'normal');

%% Plot 11: Distribution of Total Deviations
subplot(3, 4, 11);
histogram(E_total, 30, 'Normalization', 'probability', ...
    'FaceColor', [0.8, 0.4, 0.3], 'EdgeColor', 'k');
xlabel('E_{total}^{(i)}');
ylabel('Probability');
title('Distribution of Total Deviations');
grid on;

%% Plot 12: Sample Eigenfunction Trajectories
subplot(3, 4, 12);
k_plot = 1;  % First selected mode
hold on;
for m = [1, 25, 50, 75, 100]
    plot(squeeze(PHI_matrix(k_plot, :, m)), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('enPREC%03d', m-1));
end
plot(PHI_mean_matrix(k_plot, :), 'k--', 'LineWidth', 2.5, 'DisplayName', 'Mean');
xlabel('Time t');
ylabel(['\phi_{', num2str(eig_indices(k_plot)), '}(x(t))']);
title(['Eigenfunction Trajectories: Mode ', num2str(eig_indices(k_plot))]);
legend('Location', 'best', 'FontSize', 8);
grid on;

sgtitle({'\textbf{Dual Koopman Deviation Dynamics}', ...
    '$\epsilon_{t,\ell}^{(i)} = \sum_j \lambda_j^t [\phi_j(x_0^{(i)}) - \phi_j(\bar{x}_0)] \langle \zeta_j, h_\ell \rangle$'}, ...
    'Interpreter', 'latex', 'FontSize', 14);

%% ============================================================
%% STEP 9: Summary Statistics
%% ============================================================
fprintf('\n========================================\n');
fprintf('DUAL KOOPMAN DEVIATION ANALYSIS SUMMARY\n');
fprintf('========================================\n');
fprintf('Ensemble members: %d (enPREC000-099)\n', M);
fprintf('Mean trajectory: enPREC100\n');
fprintf('Time steps analyzed: %d\n', min_size);
fprintf('Spatial dimensions: %d\n', nSpace);
fprintf('Modes analyzed: %d\n\n', K);

fprintf('EIGENVALUE SPECTRUM:\n');
for k = 1:K
    fprintf('  λ_%d = %.6f + %.6fi (|λ| = %.6f)\n', ...
        eig_indices(k), real(lambda(eig_indices(k))), ...
        imag(lambda(eig_indices(k))), abs(lambda(eig_indices(k))));
end
fprintf('\n');

fprintf('EIGENFUNCTION DEVIATION ||Δφ_j^(i)||_2:\n');
for k = 1:K
    fprintf('  Mode φ_%d:\n', eig_indices(k));
    fprintf('    Mean: %.4e, Std: %.4e, Max: %.4e\n', ...
        mean(E_phi(k, :)), std(E_phi(k, :)), max(E_phi(k, :)));
end
fprintf('\n');

fprintf('NORMALIZED DEVIATION Ẽ_j^(i):\n');
for k = 1:K
    fprintf('  Mode φ_%d:\n', eig_indices(k));
    fprintf('    Mean: %.4e, Std: %.4e, Max: %.4e\n', ...
        mean(E_phi_normalized(k, :)), std(E_phi_normalized(k, :)), ...
        max(E_phi_normalized(k, :)));
end
fprintf('\n');

fprintf('TOTAL DEVIATION ENERGY:\n');
fprintf('  Mean: %.4e\n', mean(E_total));
fprintf('  Median: %.4e\n', median(E_total));
fprintf('  Std: %.4e\n', std(E_total));
fprintf('  Max: %.4e (ensemble %d = enPREC%03d)\n', ...
    max(E_total), find(E_total == max(E_total), 1)-1, ...
    find(E_total == max(E_total), 1)-1);
fprintf('\n');

fprintf('STATE-SPACE DEVIATION:\n');
fprintf('  Mean ||ε^(i)||_F: %.4e\n', mean(epsilon_norm));
fprintf('  Max ||ε^(i)||_F: %.4e (ensemble %d)\n', ...
    max(epsilon_norm), find(epsilon_norm == max(epsilon_norm), 1)-1);
fprintf('\n');

fprintf('PAIRWISE DEVIATIONS (mean across all pairs):\n');
for k = 1:K
    fprintf('  Mode φ_%d: D̄ = %.4e\n', eig_indices(k), D_pairwise_mean(k));
end
fprintf('\n');

% fprintf('CORRELATION:\n');
% fprintf('  State dev vs Eigenfunction dev: ρ = %.3f (p = %.2e)\n', rho, pval);
% fprintf('\n');

fprintf('INTERPRETATION:\n');
fprintf('  • Ẽ_j^(i): Normalized temporal deviation of mode j in ensemble i\n');
fprintf('  • E_total^(i): Combined deviation across all modes\n');
fprintf('  • γ_j^(i): Percentage contribution of mode j to total deviation\n');
fprintf('  • D_k^(i,j): Pairwise distance between ensembles i and j for mode k\n');
fprintf('========================================\n');

fprintf('\nAnalysis complete! ✓\n');