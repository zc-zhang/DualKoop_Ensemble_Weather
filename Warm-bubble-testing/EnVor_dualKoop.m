%% Ensemble data collection

% load('D:\Susuki Lab\Testing_Code\data-weather\Ensemble SCALE data Test\Ensemble-DMDs\Ensembledata20240802\data_Vfull')

% alos see load ('D:\Susuki Lab\Testing_Code\data-weather\Data250525\Dual Koopman-Ensemble-260205\Warm-bubble-testing\EnVorticityData')

% ensemble vorticity data each Vfull%3.mat: 3880*121, ensemble M= 10 ensmeble:
% p =38880, Tsnap=121, M =10;
% here we consider the transient behavior by truncated Tnsap=51; 

%% Ensemble data generated
% % Define paths and parameters
% baseDir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Dual Koopman-Ensemble-260205\Warm-bubble-testing\EnVorticityData\';
% savePathX = fullfile(baseDir, 'EnVfull_Xa.mat');
% savePathY = fullfile(baseDir, 'EnVfull_Ya.mat');
% savePath = fullfile(baseDir, 'EnVfull_all260708.mat');
% 
% M = 10;          % Number of ensembles
% p = 3880;        % Number of spatial points
% Tsnap = 51;      % Base snapshot truncation limit
% N_sub = Tsnap - 1; % 50 snapshots per ensemble for X and Y
% 
% % Preallocate the matrices: 3880 x 500
% EnVfull_Xa = zeros(p, M * N_sub);
% EnVfull_Ya = zeros(p, M * N_sub);
% 
% % Loop through each ensemble file
% for i = 0:M-1
%     % Construct the file name (Vfull000.mat to Vfull009.mat)
%     fileName = sprintf('Vfull%03d.mat', i);
%     filePath = fullfile(baseDir, fileName);
% 
%     fprintf('Processing: %s\n', fileName);
% 
%     % Load the file
%     data = load(filePath);
% 
%     % Extract the matrix safely
%     if isfield(data, 'Vfull')
%         current_matrix = data.Vfull;
%     else
%         varNames = fieldnames(data);
%         current_matrix = data.(varNames{1});
%     end
% 
%     % Extract X (1 to 50) and Y (2 to 51) snapshots for the current ensemble
%     X_ensemble = current_matrix(:, 1:N_sub);
%     Y_ensemble = current_matrix(:, 2:Tsnap);
% 
%     % Calculate column indices for insertion into the master matrices
%     startCol = (i * N_sub) + 1;
%     endCol   = (i + 1) * N_sub;
% 
%     % Insert into master matrices
%     EnVfull_Xa(:, startCol:endCol) = X_ensemble;
%     EnVfull_Ya(:, startCol:endCol) = Y_ensemble;
% end
% 
% % Save both matrices separately
% save(savePathX, 'EnVfull_Xa', '-v7.3');
% save(savePathY, 'EnVfull_Ya', '-v7.3');
% 
% fprintf('\nSuccessfully saved:\n1. %s\n2. %s\n', savePathX, savePathY);
% save(savePath, 'EnVfull_all', '-v7.3'); 
% fprintf('Successfully saved combined dataset to: %s\n', savePath);

%clc; clear all

%% Load (EnVfull_all260708.mat)
Xa_all =EnVfull_Xa;
Ya_all =EnVfull_Ya;


p=size(Xa_all,1);
M=10;
T_snap=size(Xa_all,2)/M;

N_dict = 5; % Choose your dictionary size (e.g., 150 features)
fprintf('Running kernel_ResDMD with N = %d features...\n', N_dict);

%%Run the original code signature to get the clean truncated spaces
%[G, K_star, L, PX, PY, PSI_x, PSI_y, PSI_y2, G1, A1, kernel_f] = ...
   % kernel_ResDMD(Xa_all, Ya_all, 'type', 'Gaussian', 'N', N_dict);
   [G, K_star, L, PX, PY, G1, kernel_f] = ...
    kernel_ResDMD(Xa_all, Ya_all, 'type', 'Gaussian', 'N', N_dict);

% 3. COMPUTE THE DISCRETE FEATURE COMPONENTS HERE (FIXED)
fprintf('Computing KEs, KEFs, and KMs on the truncated space...\n');

%% ========================= NOn TRunction ============
    % No truncation needed — G, K_star, PX are already N_dict-sized
    [V_coeff, Lambda_mat] = eig(K_star, G);   % G == eye(N_dict), so this ≈ eig(K_star)
    
    KEs = diag(Lambda_mat);
    KEFs = PX * V_coeff;

   % KMs     = Xa_all * pinv(KEFs).';
   KMs = Xa_all * pinv(KEFs.'); 
   %KMs = (KEFs.' \ Xa_all.').' ;   % uses QR, much faster
    Ya_pred = KMs * Lambda_mat * KEFs.';


%%Calculate and display tracking error
reconstruction_error = norm(Ya_all - real(Ya_pred), 'fro') / norm(Ya_all, 'fro');
fprintf('=======================================\n');
fprintf('Trajectory Reconstruction Error: %.4f%%\n', reconstruction_error * 100);
fprintf('=======================================\n');



%%  VISUALIZATION 1: True vs. Reconstructed Attractor Shapes
figure('Position', [100, 100, 950, 420]);

%%Left Panel: True Dynamics
subplot(1, 2, 1);
scatter3(Ya_all(1,:), Ya_all(2,:), Ya_all(3,:), 12, Ya_all(3,:), 'filled');
title('Lorenz Attractor (Ground Truth)');
xlabel('X'); ylabel('Y'); zlabel('Z');
view(45, 20); grid on; colorbar;

%%Right Panel: Reconstructed Dynamics
subplot(1, 2, 2);
scatter3(Ya_pred(1,:), Ya_pred(2,:), Ya_pred(3,:), 12, Ya_pred(3,:), 'filled');
%scatter3(Ya_all(1,:), Ya_all(2,:), Ya_all(3,:), 12, real(Ya_pred(3,:)), 'filled');
title(['Reconstructed Attractor (Dual Koopman, N_d = ' num2str(N_dict) ')']);
xlabel('X'); ylabel('Y'); zlabel('Z');
view(45, 20); grid on; colorbar;

% 6. VISUALIZATION 2: Discrete Koopman Eigenvalues Spectrum
figure;
theta = linspace(0, 2*pi, 100);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1.5); hold on;
scatter(real(KEs), imag(KEs), 45, 'r', 'filled');
title('Koopman Eigenvalue Distribution');
xlabel('Real(\lambda)'); ylabel('Imag(\lambda)');
grid on; axis equal;
legend('Unit Circle', 'Calculated KEs');

























































%% Additional RESIDUAL compute 
%% --- 1. Compute the Residuals (Your Current Block) ---
% denominators = sum(abs(V_coeff).^2, 1).'; 
% numerators = real(sum(conj(V_coeff) .* (L * V_coeff), 1)).';
% RES = sqrt(max(0, (numerators ./ denominators) - abs(KEs).^2));

RES = zeros(N_dict, 1);
for j = 1:N_dict
    v   = V_coeff(:, j);
    lam = KEs(j);
    % Full residual: ||Ug - λg||² / ||g||²
    res_num = real(v' * (L - conj(lam)*K_star - lam*K_star' + abs(lam)^2*G) * v);
    res_den = real(v' * G * v);
    RES(j)  = sqrt(max(0, res_num / res_den));
end
%[~, sort_res_idx] = sort(RES, 'ascend');  % small residual = genuine eigenvalue

%% --- 2. Sort by Residual in ASCENDING Order (Best Physics First) ---
[sorted_RES, sort_res_idx] = sort(RES, 'ascend');

%% --- 3. Evaluate Reconstruction Across Your Truncation Set ---
%%  Compute Trajectory Errors ---
N_dict_trunc_set = [10, 20, 40, 60, 90, 300]; 
num_tests = length(N_dict_trunc_set);
Ya_all_tensor = reshape(real(Ya_all), p, T_snap, M);
error_data = zeros(M, num_tests);

% Initialize cell array to store the 3D tensor for each test size
Ya_pred_N_ditcs = cell(num_tests, 1);

for i = 1:num_tests
    active_indices = sort_res_idx(1:N_dict_trunc_set(i));
    
    % Quick Reconstruction Slicing
    V_trunc = V_coeff(:, active_indices);
    L_trunc = diag(KEs(active_indices));
    Ya_pred_t = (Xa_all * pinv(PX * V_trunc).') * L_trunc * (PX * V_trunc).';
    Ya_pred_tensor = reshape(real(Ya_pred_t), 3, T_snap, M);
    
    Ya_pred_N_ditcs{i} = Ya_pred_tensor;
    % Track individual trajectory errors
    for m = 1:M
        error_data(m, i) = (norm(Ya_all_tensor(:,:,m) - Ya_pred_tensor(:,:,m), 'fro') / ...
                            norm(Ya_all_tensor(:,:,m), 'fro')) * 100;
    end
end

%% --- 2. Plot Both Side-by-Side (Ultra Simple) ---
figure('Color', 'w', 'Position', [100, 100, 1100, 500]);
x_labels = arrayfun(@(x) ['N=', num2str(x)], N_dict_trunc_set, 'UniformOutput', false);

% LEFT PANEL: Simple Standard Boxplot
%subplot(1, 2, 1);
boxplot(error_data, 'Labels', x_labels, 'Widths', 0.5);
title('Boxplot Error Distribution', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Relative Reconstruction Error (%)', 'FontSize', 11);
xlabel('Active RKHS Features (N)', 'FontSize', 11);
ylim([0, 105]);
grid on;



%% --- Define Which Loop Indices You Want to Plot ---
idx_Ndict1 = 5; % Pulls the 2nd dictionary size from your set (e.g., N = 20)
idx_Ndict2 = 6; % Pulls the 6th dictionary size from your set (e.g., N = 100)

%% --- Widescreen Plotting via Direct Cell Extraction ---
track_colors = jet(M); 

figure('Color', 'w', 'Position', [100, 100, 1200, 550]);

% =========================================================================
% SUBPLOT 1: Ground Truth
% =========================================================================
subplot(1, 3, 1);
hold on;
for m = 1:M
    % 1. Plot the trajectory line
    plot3(Ya_all_tensor(1,:,m), Ya_all_tensor(2,:,m), Ya_all_tensor(3,:,m), ...
          'Color', [track_colors(m, :), 0.4], 'LineWidth', 1.5); 
    % 2. Highlight the Initial Condition (Time Step 1)
    plot3(Ya_all_tensor(1,1,m), Ya_all_tensor(2,1,m), Ya_all_tensor(3,1,m), ...
          'o', 'MarkerFaceColor', track_colors(m, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
end

% --- ADDED: Compute and Overlay True Ensemble Mean ---
Ya_mean_true = mean(Ya_all_tensor, 3); % Size: [3 x T_snap]
plot3(Ya_mean_true(1,:), Ya_mean_true(2,:), Ya_mean_true(3,:), 'k--o', ...
      'LineWidth', 2.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k');

title('True Ensemble (Ground Truth)', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('X'); ylabel('Y'); zlabel('Z');
view(45, 20); grid on;
hold off;

% =========================================================================
% SUBPLOT 2: First Selected Dictionary Truncation
% =========================================================================
subplot(1, 3, 2);
hold on;
for m = 1:M
    % 1. Plot the reconstructed trajectory line
    plot3(Ya_pred_N_ditcs{idx_Ndict1}(1,:,m), Ya_pred_N_ditcs{idx_Ndict1}(2,:,m), Ya_pred_N_ditcs{idx_Ndict1}(3,:,m), ...
          'Color', [track_colors(m, :), 0.4], 'LineWidth', 1.5); 
    % 2. Highlight the Initial Condition (Time Step 1)
    plot3(Ya_pred_N_ditcs{idx_Ndict1}(1,1,m), Ya_pred_N_ditcs{idx_Ndict1}(2,1,m), Ya_pred_N_ditcs{idx_Ndict1}(3,1,m), ...
          'o', 'MarkerFaceColor', track_colors(m, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
end

% --- ADDED: Compute and Overlay Reconstructed Mean (Truncated) ---
Ya_mean_pred1 = mean(Ya_pred_N_ditcs{idx_Ndict1}, 3);
plot3(Ya_mean_pred1(1,:), Ya_mean_pred1(2,:), Ya_mean_pred1(3,:), 'k--o', ...
      'LineWidth', 2.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k');
% -----
title(['Reconstructed (N = ' num2str(N_dict_trunc_set(idx_Ndict1)) ')'], 'FontSize', 11, 'FontWeight', 'bold');
xlabel('X'); ylabel('Y'); zlabel('Z');
view(45, 20); grid on;
hold off;

% =========================================================================
% SUBPLOT 3: Second Selected Dictionary Truncation
% =========================================================================
subplot(1, 3, 3);
hold on;
for m = 1:M
    % 1. Plot the reconstructed trajectory line
    plot3(Ya_pred_N_ditcs{idx_Ndict2}(1,:,m), Ya_pred_N_ditcs{idx_Ndict2}(2,:,m), Ya_pred_N_ditcs{idx_Ndict2}(3,:,m), ...
          'Color', [track_colors(m, :), 0.4], 'LineWidth', 1.5); 
    % 2. Highlight the Initial Condition (Time Step 1)
    plot3(Ya_pred_N_ditcs{idx_Ndict2}(1,1,m), Ya_pred_N_ditcs{idx_Ndict2}(2,1,m), Ya_pred_N_ditcs{idx_Ndict2}(3,1,m), ...
          'o', 'MarkerFaceColor', track_colors(m, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
end

% --- ADDED: Compute and Overlay Reconstructed Mean (Full Dictionary) ---
Ya_mean_pred2 = mean(Ya_pred_N_ditcs{idx_Ndict2}, 3);
plot3(Ya_mean_pred2(1,:), Ya_mean_pred2(2,:), Ya_mean_pred2(3,:), 'k--o', ...
      'LineWidth', 2.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k');

title(['Reconstructed (N = ' num2str(N_dict_trunc_set(idx_Ndict2)) ')'], 'FontSize', 11, 'FontWeight', 'bold');
xlabel('X'); ylabel('Y'); zlabel('Z');
view(45, 20); grid on;
hold off;
