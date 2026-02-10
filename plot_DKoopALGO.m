% plot_DKoopALGO.m    % shortenhand
% plot_dualKoop_Algom.m  % full name

% plot_DKoopALGO_corrected.m
data_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Ensemble_PREC';
M = 101;
Xa_ens = cell(M,1);
Ya_ens = cell(M,1);

for m = 0:M-1
    fname = sprintf('enPREC%03d.mat', m);
    S = load(fullfile(data_dir, fname));
    Xraw = S.vectorized_PREC;    % 20860 × 72
    
    % DON'T TRANSPOSE! Keep same format as single trajectory
    Xa_ens{m+1} = Xraw(:,2:end-1);  % 20860 × 72 
    Ya_ens{m+1} = Xraw(:,3:end);    % 20860 × 72 
end


%% Choose reference ensemble
Xa_ref = Xa_ens{101};
Ya_ref = Ya_ens{101};

%% Build operator ONCE using reference ensemble
%fprintf('Building dual Koopman operator on reference ensemble...\n');
%% Build ONCE
[G, K_star, L, lambda, W, C_ref, PSI_ref, ~, ~, ~, UU_ref] = dualKoop_Algom( ...
    Xa_ref, Ya_ref, ...
    'type', 'Linear', ...
    'N', 300, 'cut_off', ...
    1e-10 ...
);

%% Evaluate efficiently
parfor m = 1:M
    [~,~,~,~,~,~,PSI_ens{m}] = dualKoop_Algom( ...
        Xa_ref, Ya_ref, ...
        'Xb', Xa_ens{m}, ...
        'UU', UU_ref ...      % ← Reuse! No re-eigendecomposition
    );

end

%% Extract eigenvectors for projection (if available from your function)
% You need to store the eigenvectors UU from the reference computation
% This assumes dualKoop_Algom can return them or you modify it

% CRITICAL: Use the SAME kernel parameters for all ensembles
% Store kernel bandwidth from reference
d_ref = mean(vecnorm(Xa_ref - mean(Xa_ref, 2)));

%% Evaluate dual eigenfunctions on each ensemble
PSI_ens = cell(M, 1);
PSI_ens{101} = PSI_ref;  % Store reference

fprintf('Evaluating eigenfunctions on %d ensembles...\n', M);
parfor m = 1:M
    if m == 101
        continue;  % Already have reference
    end
    
    % Option 1: If dualKoop_Algom supports external data evaluation
    [~, ~, ~, ~, ~, PSI_ens{m}] = dualKoop_Algom( ...
        Xa_ref, Ya_ref, ...  % Use REFERENCE for operator construction
        'type', 'Gaussian', ...
        'N', 300, ...
        'cut_off', 1e-10, ...
        'Xb', Xa_ens{m} ...  % Evaluate on this ensemble
    );
end

%% Handle size differences
min_size = min(cellfun(@(x) size(x, 1), PSI_ens));
PSI_aligned = cellfun(@(P) P(1:min_size, :), PSI_ens, 'UniformOutput', false);

%% Select eigenfunction
eig_id = 2;

%% Plot eigenvalue spectrum first
figure;
semilogy(1:length(lambda), abs(lambda), 'o-', 'LineWidth', 1.5);
xlabel('Eigenfunction index');
ylabel('|Eigenvalue|');
title('Koopman Eigenvalue Spectrum');
grid on;

%% Check which are significant
fprintf('Top 10 eigenvalues:\n');
for k = 1:min(10, length(lambda))
    fprintf('  eig_id = %d: λ = %.6f (|λ| = %.6f)\n', k, lambda(k), abs(lambda(k)));
end

%% Stack all eigenfunctions into matrix
phi_all = zeros(min_size, M-1);
parfor m = 1:M-1
    phi_m = PSI_aligned{m}(:, eig_id);
    %phi_all(:, m) = phi_m / norm(phi_m);  % Normalize while stacking
    phi_all(:, m) = phi_m ;  % Normalize while stacking
end

%% Compute mean
phi_mean = mean(phi_all, 2);  % min_size × 1
%phi_mean = phi_mean / norm(phi_mean);


phi_ref = phi_all(:, 100);  % Or use mean, your choice

%% Compute errors
err_to_ref = zeros(M-1, 1);
err_to_mean = zeros(M-1, 1);

parfor m = 1:M-1
    phi_m = phi_all(:, m);
    err_to_ref(m) = norm(phi_m - phi_ref) / norm(phi_ref);
    err_to_mean(m) = norm(phi_m - phi_mean) / norm(phi_mean);
end

%% Plot
figure;
semilogy(err_to_mean, 'o-', 'LineWidth', 1.5);
xlabel('Ensemble index m');
ylabel('Relative discrepancy from mean');
title(['Ensemble consistency of dual KEF \phi_', num2str(eig_id)]);
grid on;

%% Visualization
figure('Position', [100, 100, 1200, 400]);

% Plot 1: Discrepancy to reference
subplot(1, 3, 1);
semilogy(err_to_ref, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('Ensemble index m');
ylabel('Relative discrepancy');
title(['Discrepancy to Reference (enPREC100), \phi_', num2str(eig_id)]);
grid on;
ylim([0, max(err_to_ref)*1.1]);

% Plot 2: Discrepancy to ensemble mean
subplot(1, 3, 2);
plot(0:M-1, err_to_mean, 's-', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('Ensemble index m');
ylabel('Relative discrepancy');
title(['Discrepancy to Ensemble Mean, \phi_', num2str(eig_id)]);
grid on;
ylim([0, max(err_to_mean)*1.1]);

% Plot 3: Eigenfunction samples
%% Plot: Sample eigenfunctions over time
%% Align signs to ensemble mean BEFORE plotting
for m = 1:M
    phi_m = phi_all(:, m);
    if dot(phi_m, phi_mean) < 0
        phi_all(:, m) = -phi_m;  % Flip sign
    end
end

%% NOW plot with aligned signs
subplot(1, 3, 3);
hold on;

for m = [1, 25, 50, 75, 101]
    phi_m = phi_all(:, m);  % Use phi_all (already aligned)
    plot(real(phi_m), 'DisplayName', sprintf('enPREC%03d', m-1), 'LineWidth', 1.2);
end

plot(real(phi_mean), 'k--', 'LineWidth', 2, 'DisplayName', 'Ensemble Mean');
legend('Location', 'best');
xlabel('Time index');
ylabel(['\phi_', num2str(eig_id), '(x(t))']);
title(['Dual Koopman Eigenfunction \phi_', num2str(eig_id)]);
grid on;
hold off;

%% Statistics
fprintf('\n=== Ensemble Consistency Analysis ===\n');
fprintf('Eigenfunction ID: %d\n', eig_id);
fprintf('Number of ensembles: %d\n', M);
fprintf('\nDiscrepancy to reference (enPREC100):\n');
fprintf('  Mean: %.4e\n', mean(err_to_ref));
fprintf('  Std:  %.4e\n', std(err_to_ref));
fprintf('  Max:  %.4e (ensemble %d)\n', max(err_to_ref), find(err_to_ref == max(err_to_ref))-1);
fprintf('\nDiscrepancy to ensemble mean:\n');
fprintf('  Mean: %.4e\n', mean(err_to_mean));
fprintf('  Std:  %.4e\n', std(err_to_mean));
fprintf('  Max:  %.4e (ensemble %d)\n', max(err_to_mean), find(err_to_mean == max(err_to_mean))-1);

%% Additional diagnostic: Check if discrepancy correlates with data deviation
data_dev = zeros(M-1, 1);
for m = 1:M-1
    data_dev(m) = norm(Xa_ens{m} - Xa_ref, 'fro') / norm(Xa_ref, 'fro');
end

figure;
scatter(data_dev, err_to_ref, 50, 'filled');
xlabel('Data deviation from reference (Frobenius norm)');
ylabel('Eigenfunction discrepancy');
title('Correlation: Data vs Eigenfunction Deviation');
grid on;
[rho, pval] = corr(data_dev, err_to_ref);
text(0.05, 0.95, sprintf('\\rho = %.3f (p = %.3e)', rho, pval), ...
    'Units', 'normalized', 'FontSize', 12);



%%
%% Compute pairwise discrepancies
D_pairwise = zeros(M-1, M-1);

for i = 1:M-1
    for j = 1:M-1
        if i == j
            D_pairwise(i,j) = 0;
        else
            phi_i = phi_all(:, i);
            phi_j = phi_all(:, j);
            D_pairwise(i,j) = norm(phi_i - phi_j) / norm(phi_j);
        end
    end
end

%% Visualize pairwise matrix
figure;
imagesc(D_pairwise);
colorbar;
xlabel('Ensemble j');
ylabel('Ensemble i');
title(['Pairwise Discrepancy Matrix for \phi_', num2str(eig_id)]);
axis equal tight;




%%%%%%%% New check plot 

% %% Extract eigenfunction for selected ensembles
% eig_id = 3;  % Try different values: 1, 2, 3, ...
% 
% % Build matrix of all eigenfunctions
% phi_all = zeros(min_size, M);
% for m = 1:M-1
%     phi_m = PSI_aligned{m}(:, eig_id);
%     phi_all(:, m) = phi_m / norm(phi_m);  % Normalize
% end
% 
% % Compute mean
% phi_mean = mean(phi_all, 2);
% phi_mean = phi_mean / norm(phi_mean);
% 
% %% ALIGN SIGNS relative to mean
% for m = 1:M-1
%     if dot(phi_all(:, m), phi_mean) < 0
%         phi_all(:, m) = -phi_all(:, m);
%     end
% end
% 
% % Recompute mean after alignment
% phi_mean = mean(phi_all, 2);
% phi_mean = phi_mean / norm(phi_mean);
% 
% %% Plot
% figure('Position', [100, 100, 1400, 500]);
% 
% % Subplot 1: Overlay
% subplot(1, 3, 1);
% hold on;
% for m = [1, 25, 50, 75, 101]
%     plot(phi_all(:, m), 'DisplayName', sprintf('enPREC%03d', m-1), 'LineWidth', 1.2);
% end
% plot(phi_mean, 'k--', 'LineWidth', 2.5, 'DisplayName', 'Mean');
% legend('Location', 'best');
% xlabel('Time index');
% ylabel(['\phi_', num2str(eig_id), '(x(t))']);
% title(['Eigenfunction \phi_', num2str(eig_id), ' (λ=', num2str(lambda(eig_id), '%.4f'), ')']);
% grid on;
% 
% % Subplot 2: Standard deviation band
% subplot(1, 3, 2);
% phi_std = std(phi_all, 0, 2);  % Std across ensembles
% t = 1:min_size;
% fill([t, fliplr(t)], [phi_mean' + phi_std', fliplr(phi_mean' - phi_std')], ...
%     'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '±1 std');
% hold on;
% plot(phi_mean, 'k-', 'LineWidth', 2, 'DisplayName', 'Mean');
% xlabel('Time index');
% ylabel(['\phi_', num2str(eig_id), '(x(t))']);
% title('Mean ± Standard Deviation');
% legend;
% grid on;
% 
% % Subplot 3: Discrepancies
% subplot(1, 3, 3);
% err_to_mean = zeros(M, 1);
% for m = 1:M-1
%     err_to_mean(m) = norm(phi_all(:, m) - phi_mean) / norm(phi_mean);
% end
% plot(0:M-1, err_to_mean, 'o-', 'LineWidth', 1.5);
% xlabel('Ensemble index');
% ylabel('Relative error to mean');
% title('Ensemble Consistency');
% grid on;
% 
% fprintf('\n=== Eigenfunction %d Analysis ===\n', eig_id);
% fprintf('Eigenvalue: λ = %.6f (|λ| = %.6f)\n', lambda(eig_id), abs(lambda(eig_id)));
% fprintf('Mean discrepancy: %.4e\n', mean(err_to_mean));
% fprintf('Std discrepancy:  %.4e\n', std(err_to_mean));
% fprintf('Max discrepancy:  %.4e (ensemble %d)\n', max(err_to_mean), find(err_to_mean == max(err_to_mean))-1);