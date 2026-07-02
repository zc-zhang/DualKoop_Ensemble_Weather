%% 1. LOAD YOUR ENSEMBLE DATA
% load('LZ_ensemble_data_dim3.mat'); % Contains Xa (3 x M) and Ya (3 x M)
% [dim, M] = size(Xa);


%% 1-1 generated the LZ ensmeble data 

%% ====================================================================
%  COMPLETE INTEGRATED DUAL-KOOPMAN ATTRACTOR RECONSTRUCTION
%  REQUIRES: kernel_ResDMD.m (unmodified) on the MATLAB path.
% ====================================================================
clear; close all; rng(1);

%% -------------------- STEP 1: Lorenz Ensemble Generation --------------------
SIGMA = 10; RHO = 28; BETA = 8/3;   % classic chaotic Lorenz parameters
ODEFUN = @(t,y) [SIGMA*(y(2)-y(1));
                 y(1).*(RHO-y(3))-y(2);
                 y(1).*y(2)-BETA*y(3)];
             
M       = 50;     % number of ensemble members (initial conditions)
delta_t = 0.05;   % sampling interval for X->Y snapshot pairs
T_snap  = 65;     % number of snapshot pairs PER ensemble member
MT_snap = M*T_snap;
T_burn  = 5;    % burn-in time so each IC is already near/on the attractor
options = odeset('RelTol',1e-10,'AbsTol',1e-11);

%X0_raw  = 1*(rand(3,M)-0.5) + [0;0;10];   % crude spread of raw ICs
% X0_raw  = 1*(rand(3,M)-0.01) + [0;0;10];   % crude spread of raw ICs
% Xa_ens  = cell(M,1);  Ya_ens = cell(M,1);  X0_all = zeros(3,M);
% 
% fprintf('Generating %d ensemble Lorenz trajectories (ICs vary only)...\n', M);
% for m = 1:M
%     ode_use = ODEFUN;
%     [~,Yb] = ode45(ode_use, [0 T_burn], X0_raw(:,m), options);
%     x0 = Yb(end,:)';  X0_all(:,m) = x0;
%     t_grid = 0:delta_t:(T_snap*delta_t);
%     [~,Ytraj] = ode45(ode_use, t_grid, x0, options);
%     Ytraj = Ytraj';
%     Xa_ens{m} = Ytraj(:,1:end-1);
%     Ya_ens{m} = Ytraj(:,2:end);
% end

%% CRITICAL FIX: Initialize all members in a tight ball on ONE WING center
% UNIQUE LOCALIZED BATCH: Pick a point on the outer edge of the attractor
% [x, y, z] = [-10, -15, 25] is a great launchpad that forces a full two-wing split later
base_IC = [-10.0; -15.0; 25.0]; 

% A very small perturbation (scaled by 0.01) creates an ultra-tight localized cluster
X0_raw = 0.015 * (rand(3, M) - 0.5) + base_IC;

Xa_ens  = cell(M,1);  Ya_ens = cell(M,1);  X0_all = zeros(3,M);

for m = 1:M
    % Brief burn-in to settle the tight batch onto the attractor manifold smoothly
    [~,Y_burn] = ode45(ODEFUN, [0 T_burn], X0_raw(:,m), options);
    X0 = Y_burn(end,:).';
    X0_all(:,m) = X0;
    
    % Core simulation tracking
    t_span = 0:delta_t:(T_snap*delta_t);
    [~,Y_orbit] = ode45(ODEFUN, t_span, X0, options);
    
    Xa_ens{m} = Y_orbit(1:end-1,:).'; 
    Ya_ens{m} = Y_orbit(2:end,:).';   
end

% Pool ALL members for dual dictionary mapping
Xa_all = cell2mat(Xa_ens');
Ya_all = cell2mat(Ya_ens');
[dim, Total_Snapshots] = size(Xa_all);
fprintf('Done. Total snapshots pooled: %d\n', Total_Snapshots);


%% --- SAVE POOLED DATA TO FILE ---
save('LZ_ensemble_raw_data01.mat', 'Xa_all', 'Ya_all');
fprintf('Pooled ensemble data successfully saved to LZ_ensemble_raw_data01.mat!\n');

% % %% 2. RUN COLBROOK'S ORIGINAL SCRIPT WITH TRUNCATION SIZE N
N_dict = 100; % Choose your dictionary size (e.g., 150 features)
fprintf('Running kernel_ResDMD with N = %d features...\n', N_dict);

%%Run the original code signature to get the clean truncated spaces
[G, K_star, L, PX, PY, PSI_x, PSI_y, PSI_y2, G1, A1, kernel_f] = ...
    kernel_ResDMD(Xa_all, Ya_all, 'type', 'Gaussian', 'N', N_dict);

% 3. COMPUTE THE DISCRETE FEATURE COMPONENTS HERE (FIXED)
fprintf('Computing KEs, KEFs, and KMs on the truncated space...\n');

%% ========================= NOn TRunction ============
    % No truncation needed — G, K_star, PX are already N_dict-sized
    [V_coeff, Lambda_mat] = eig(K_star, G);   % G == eye(N_dict), so this ≈ eig(K_star)
    
    KEs = diag(Lambda_mat);
    KEFs = PX * V_coeff;

    KMs     = Xa_all * pinv(KEFs).';
    Ya_pred = KMs * Lambda_mat * KEFs.';


%%Calculate and display tracking error
reconstruction_error = norm(Ya_all - Ya_pred, 'fro') / norm(Ya_all, 'fro');
fprintf('=======================================\n');
fprintf('Trajectory Reconstruction Error: %.4f%%\n', reconstruction_error * 100);
fprintf('=======================================\n');


%% Additional RESIDUAL compute 
%% --- 1. Compute the Residuals (Your Current Block) ---
denominators = sum(abs(V_coeff).^2, 1).'; 
numerators = real(sum(conj(V_coeff) .* (L * V_coeff), 1)).';
RES = sqrt(max(0, (numerators ./ denominators) - abs(KEs).^2));

%% --- 2. Sort by Residual in ASCENDING Order (Best Physics First) ---
[sorted_RES, sort_residual_idx] = sort(RES, 'ascend');

%% --- 3. Evaluate Reconstruction Across Your Truncation Set ---
%%  Compute Trajectory Errors ---
N_dict_trunc_set = [10, 20, 40, 60, 80, 100]; 
num_tests = length(N_dict_trunc_set);
Ya_all_tensor = reshape(real(Ya_all), 3, T_snap, M);
error_data = zeros(M, num_tests);

% Initialize cell array to store the 3D tensor for each test size
Ya_pred_N_ditcs = cell(num_tests, 1);

for i = 1:num_tests
    active_indices = sort_residual_idx(1:N_dict_trunc_set(i));
    
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
idx_Ndict1 = 2; % Pulls the 2nd dictionary size from your set (e.g., N = 20)
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
title(['Reconstructed (N = ' num2str(N_dict_trunc_set(idx_Ndict2)) ')'], 'FontSize', 11, 'FontWeight', 'bold');
xlabel('X'); ylabel('Y'); zlabel('Z');
view(45, 20); grid on;
hold off;
%%   USE TRUNCATION IN FACT we need not use this 
% %%Extract the explicitly truncated blocks matching N_dict
% G_trunc = G(1:N_dict, 1:N_dict);
% K_trunc = K_star(1:N_dict, 1:N_dict);
% 
% %%Solve the generalized eigenvalue problem over the dictionary space
% [V_coeff, Lambda_mat] = eig(K_trunc, G_trunc);
% 
% %%1. Koopman Eigenvalues (KEs)
% KEs = diag(Lambda_mat);
% 
% %%2. Koopman Eigenfunctions (KEFs) evaluated at the training snapshots
% %%Since PX = G1 * UU, extracting the first N_dict columns of PX gives exactly G1 * UU_trunc!
% PX_trunc = PX(:, 1:N_dict); 
% KEFs = PX_trunc * V_coeff; 
% 
% %%3. Koopman Modes (KMs)
% %%To get KMs = Xa * UU_trunc * V_coeff, we use the identity UU = G1 \ PX
% % KMs = Xa_all * (G1 \ PX_trunc) * V_coeff;
% % 
% % Ya_pred = KMs * (Lambda_mat * pinv(KEFs) * G1);
% 
% KMs    = Xa_all * pinv(KEFs).';     % state expressed in eigenfunction coords
% Ya_pred = KMs * Lambda_mat * KEFs.'; % consistent one-step forecast


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


%% --- Standalone Figure: True Lorenz Dynamics with Ensemble Marks ---
figure('Color', 'w', 'Position', [100, 100, 650, 550]);
hold on;

% 1. Define a vibrant color palette so individual tracks are distinguishable
colors = jet(M); 

% 2. Loop through each trajectory to plot its continuous path
for m = 1:M
    % Extract the 3D coordinates for the current trajectory
    X_track = Ya_ens{m}(1, :);
    Y_track = Ya_ens{m}(2, :);
    Z_track = Ya_ens{m}(3, :);
    
    % Draw the continuous trajectory line
    plot3(X_track, Y_track, Z_track, 'Color', colors(m, :), 'LineWidth', 1.2);
    
    % 3. Extract and Highlight the Initial Condition (First Step t = 0)
    plot3(X_track(1), Y_track(1), Z_track(1), ...
          'o', 'MarkerSize', 8, ...
          'MarkerFaceColor', colors(m, :), ...
          'MarkerEdgeColor', [0.2 0.2 0.2], ...
          'LineWidth', 1.5);
end

% Formatting the Standalone 3D Physical Space
title('Lorenz Attractor with Ensembles (Ground Truth)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('X (Axis 1)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Y (Axis 2)', 'FontSize', 11, 'FontWeight', 'bold');
zlabel('Z (Axis 3)', 'FontSize', 11, 'FontWeight', 'bold');

view(45, 20); 
grid on;
box on;
hold off;


%% --- Figure: Reconstructed Koopman Dynamics ---
figure('Color', 'w', 'Position', [150, 150, 650, 550]);
hold on;

% 1. Define the identical color palette so trajectories match the true figure perfectly
colors = jet(M); 

% 2. Reshape Ya_pred from [3 x 2000] to a [3 x steps_per_traj x M] tensor for sequential plotting
% This extracts the reconstructed paths without needing cell wrappers
steps_per_traj = size(Ya_ens{1}, 2); % Should be 40 steps
Ya_pred_tensor = reshape(real(Ya_pred), 3, steps_per_traj, M);

% 3. Loop through and plot each reconstructed trajectory
for m = 1:M
    X_pred = Ya_pred_tensor(1, :, m);
    Y_pred = Ya_pred_tensor(2, :, m);
    Z_pred = Ya_pred_tensor(3, :, m);
    
    % Draw the continuous reconstructed trajectory line
    plot3(X_pred, Y_pred, Z_pred, 'Color', colors(m, :), 'LineWidth', 1.2);
    
    % Highlight the Reconstructed Initial Condition Footprint (Step t = 0)
    plot3(X_pred(1), Y_pred(1), Z_pred(1), ...
          'o', 'MarkerSize', 8, ...
          'MarkerFaceColor', colors(m, :), ...
          'MarkerEdgeColor', [0.2 0.2 0.2], ...
          'LineWidth', 1.5);
end

% Formatting the Standalone 3D Space
title(['Reconstructed Attractor (Dual Koopman, N_d = ' num2str(N_dict) ')'], 'FontSize', 13, 'FontWeight', 'bold');
xlabel('X (Axis 1)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Y (Axis 2)', 'FontSize', 11, 'FontWeight', 'bold');
zlabel('Z (Axis 3)', 'FontSize', 11, 'FontWeight', 'bold');

view(45, 20); 
grid on;
box on;
hold off;




%% SUBPLOT 

%% --- Combined Layout: Subplots with LaTeX Rendering ---
figure('Color', 'w', 'Position', [100, 100, 1200, 550]);
colors = jet(M); 
steps_per_traj = size(Ya_ens{1}, 2); 

% =========================================================================
% LEFT PANEL: True Dynamics
% =========================================================================
subplot(1, 2, 1);
hold on;

for m = 1:M
    X_track = Ya_ens{m}(1, :);
    Y_track = Ya_ens{m}(2, :);
    Z_track = Ya_ens{m}(3, :);
    
    plot3(X_track, Y_track, Z_track, 'Color', colors(m, :), 'LineWidth', 1.2);
    plot3(X_track(1), Y_track(1), Z_track(1), ...
          'o', 'MarkerSize', 8, ...
          'MarkerFaceColor', colors(m, :), ...
          'MarkerEdgeColor', [0.2 0.2 0.2], ...
          'LineWidth', 1.5);
end

% LaTeX Formatting for Left Title and Axes
title('$$\text{Lorenz Attractor with Ensembles (Ground Truth Version: } \mathbf{Y}_a\text{)}$$', ...
      'Interpreter', 'latex', 'FontSize', 13);
xlabel('$$X$$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$$Y$$', 'Interpreter', 'latex', 'FontSize', 12);
zlabel('$$Z$$', 'Interpreter', 'latex', 'FontSize', 12);

view(45, 20); grid on; box on;
hold off;

% =========================================================================
% RIGHT PANEL: Reconstructed Dynamics
% =========================================================================
subplot(1, 2, 2);
hold on;

Ya_pred_tensor = reshape(real(Ya_pred), 3, steps_per_traj, M);

for m = 1:M
    X_pred = Ya_pred_tensor(1, :, m);
    Y_pred = Ya_pred_tensor(2, :, m);
    Z_pred = Ya_pred_tensor(3, :, m);
    
    plot3(X_pred, Y_pred, Z_pred, 'Color', colors(m, :), 'LineWidth', 1.2);
    plot3(X_pred(1), Y_pred(1), Z_pred(1), ...
          'o', 'MarkerSize', 8, ...
          'MarkerFaceColor', colors(m, :), ...
          'MarkerEdgeColor', [0.2 0.2 0.2], ...
          'LineWidth', 1.5);
end

% Dynamic LaTeX Formatting for Right Title containing N_dict
title(['$$\text{Reconstructed Attractor (Dual Koopman, } N_d = ' num2str(N_dict) '\text{)}$$'], ...
      'Interpreter', 'latex', 'FontSize', 13);
xlabel('$$X$$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$$Y$$', 'Interpreter', 'latex', 'FontSize', 12);
zlabel('$$Z$$', 'Interpreter', 'latex', 'FontSize', 12);

view(45, 20); grid on; box on;
hold off;


%% -------------------- STEP 5.6: VISUALIZE KOOPMAN MODES --------------------
fprintf('Plotting Koopman Modes (KMs) structure...\n');

% Compute the total "energy" or magnitude of each mode across all 3 physical states
% KMs is size (3 x N_dict)
km_magnitudes = sqrt(sum(abs(KMs).^2, 1)); % Size: 1 x N_dict

% ------------------ PANEL 1: MODE DOMINANCE (BAR CHART) ------------------
figure('Position', [100, 100, 1000, 420]);

subplot(1, 2, 1);
bar(1:N_dict, km_magnitudes, 'FaceColor', [0.2 0.6 0.8]);
set(gca, 'YScale', 'log'); % Often best in log scale as a few modes dominate
title('Koopman Mode Dominance (Energy Spectrum)');
xlabel('Mode Index (j)'); ylabel('Mode Magnitude ||v_j|| (Log Scale)');
grid on;

% ------------------ PANEL 2: SPATIAL VECTOR STRUCTURE ------------------
% Let's plot the top 5 most dominant modes as 3D vectors
[~, dominant_indices] = sort(km_magnitudes, 'descend');
top_modes = dominant_indices(1:10);

subplot(1, 2, 2);
hold on;
colors = lines(10);

for i = 1:10
    idx = top_modes(i);
    % Extract real parts of the 3D components for the vector arrow
    vx = real(KMs(1, idx));
    vy = real(KMs(2, idx));
    vz = real(KMs(3, idx));
    
    % Draw a 3D arrow starting from the origin
    quiver3(0, 0, 0, vx, vy, vz, 'LineWidth', 2.5, 'Color', colors(i,:), 'MaxHeadSize', 0.5);
end

title('Top 10 Dominant Koopman Modes in 3D Physical Space');
xlabel('X Direction contribution'); ylabel('Y Direction contribution'); zlabel('Z Direction contribution');
grid on; axis equal; view(45, 20);
legend(cellstr(num2str(top_modes', 'Mode %d')), 'Location', 'best');


%% -------------------- STEP 5.5: VISUALIZE KEFs (MAGNITUDE & PHASE) --------------------
fprintf('Plotting Koopman Eigenfunction Magnitude and Phase...\n');

% Select which eigenfunction mode you want to look at (e.g., Mode 1, Mode 2, etc.)
mode_to_plot = 1; 

% Extract the complex values for this specific mode across all 2400 snapshots
kef_mode = KEFs(:, mode_to_plot);

% Compute Magnitude and Phase (Angle in radians)
kef_real = real(kef_mode);
kef_phase     = angle(kef_mode); % Outputs values between -pi and +pi

% Create a split figure for side-by-side comparison
figure('Position', [100, 100, 1100, 450]);

% Left Panel: Eigenfunction Magnitude Mapping
subplot(1, 2, 1);
scatter3(Ya_all(1,:), Ya_all(2,:), Ya_all(3,:), 15, kef_real, 'filled');
title(['KEF Mode ' num2str(mode_to_plot) ' Real Part Re(\varphi)']);
xlabel('X'); ylabel('Y'); zlabel('Z');
colormap(subplot(1,2,1),jet); % 'hot' or 'magma' works great for magnitude intensities
colorbar;
view(45, 20); grid on; axis tight;

% Right Panel: Eigenfunction Phase Angle Mapping
subplot(1, 2, 2);
scatter3(Ya_all(1,:), Ya_all(2,:), Ya_all(3,:), 15, kef_phase, 'filled');
title(['KEF Mode ' num2str(mode_to_plot) ' Phase Angle \angle\varphi (rad)']);
xlabel('X'); ylabel('Y'); zlabel('Z');
colormap(subplot(1,2,2), hsv); % 'hsv' is perfect for phase because it is cyclic (-pi matches +pi)
colorbar;
view(45, 20); grid on; axis tight;






% Assumed variables from your workspace: 
% KEFs (2400 x 50), M = 10 trajectories, steps_per_traj = 240

%M = 10;
steps_per_traj = T_snap;
N_modes_to_check = 3; % Look at the first few non-trivial eigenfunctions

% 1. Extract the initial condition evaluations: < phi_j, kappa_{x_0^i} >
% Size will be [M x N_dict]
init_indices = 1:steps_per_traj:MT_snap;
KEF_initials = KEFs(init_indices, :); 

% 2. Compute your Mean Centered Drift formulation: < phi_j, kappa - bar_kappa >
mean_kappa = mean(KEF_initials, 1); % 1 x N_dict
centered_drift = KEF_initials - mean_kappa; % M x N_dict

% 3. Plot the results to see if they cleanly classify the ensemble index
figure('Position', [100, 100, 950, 450]);

% Left Panel: Raw values showing how each initialization gets a unique coordinate
subplot(1, 2, 1);
imagesc(real(KEF_initials(:, 1:N_modes_to_check))');
colorbar; colormap(hsv);
set(gca, 'YTick', 1:N_modes_to_check, 'YTickLabel', cellstr(num2str((1:N_modes_to_check)', 'Mode %d')));
title('\langle\varphi_j, \kappa_{x_0^i}\rangle (Raw Init Coordinates)');
xlabel('Ensemble Trajectory Index (i)'); ylabel('Eigenfunction');

% Right Panel: Centered drift showing how they diverge from the ensemble average
subplot(1, 2, 2);
imagesc(real(centered_drift(:, 1:N_modes_to_check))');
colorbar; colormap(hsv);
set(gca, 'YTick', 1:N_modes_to_check, 'YTickLabel', cellstr(num2str((1:N_modes_to_check)', 'Mode %d')));
title('\langle\varphi_j, \kappa_{x_0^i} - \bar{\kappa}_{x_0}\rangle (Centered Drift)');
xlabel('Ensemble Trajectory Index (i)'); ylabel('Eigenfunction');



% Assumed workspace values:
% Ya_all (3 x 2400), KEFs (2400 x 50), M = 10, steps_per_traj = 240
%M = 10;
steps_per_traj = T_snap;
mode_idx = 2; % Choose which eigenfunction channel to analyze

% 1. Extract and center the initial condition features
init_indices = 1:steps_per_traj:MT_snap;
KEF_initials = KEFs(init_indices, mode_idx); % M x 1 complex vector
mean_init = mean(KEF_initials);
centered_drift = KEF_initials - mean_init;   % M x 1 complex vector

% 2. Broadcast the single initialization value across its entire trajectory track
% This creates a 2000 x 1 vector where every point in track i shares the same value
spatial_drift_vector = zeros(MT_snap, 1);
for i = 1:M
    row_range = ((i-1)*steps_per_traj + 1) : (i*steps_per_traj);
    spatial_drift_vector(row_range) = centered_drift(i);
end

% 3. Extract spatial properties for 3D visualization
drift_magnitude = abs(spatial_drift_vector);
drift_angle     = angle(spatial_drift_vector); % Phase angle in radians (-pi to +pi)

% 4. Plot side-by-side 3D mappings
figure('Position', [50, 100, 1200, 500]);

% Left Subplot: Magnitude of Initialization Drift
subplot(1, 2, 1);
scatter3(Ya_all(1,:), Ya_all(2,:), Ya_all(3,:), 15, drift_magnitude, 'filled');
title(['Ensemble Init Drift Magnitude |\langle\varphi_{' num2str(mode_idx) '}, \kappa_{x_0^i} - \bar{\kappa}_{x_0}\rangle|']);
xlabel('X'); ylabel('Y'); zlabel('Z');
colormap(subplot(1,2,1), parula); colorbar;
view(45, 20); grid on; axis tight;

% Right Subplot: Phase Angle of Initialization Drift
subplot(1, 2, 2);
scatter3(Ya_all(1,:), Ya_all(2,:), Ya_all(3,:), 15, drift_angle, 'filled');
title(['Ensemble Init Drift Phase Angle \angle\langle\varphi_{' num2str(mode_idx) '}, \kappa_{x_0^i} - \bar{\kappa}_{x_0}\rangle']);
xlabel('X'); ylabel('Y'); zlabel('Z');
colormap(subplot(1,2,2), hsv); colorbar; % Cyclic HSV is ideal for angles
view(45, 20); grid on; axis tight;


%% Let J be the number of dominant Koopman modes you want to check (e.g., first 10 modes)
J = 30; 
KEF_initials_all = KEFs(init_indices, 1:J); % Size: M x J

% Center the matrix by subtracting the column means
mean_inits = mean(KEF_initials_all, 1);
centered_matrix = KEF_initials_all - mean_inits; % Size: M x J

% Compute the Anomaly Score (L2 norm squared of the drift per trajectory)
anomaly_scores = sum(abs(centered_matrix).^2, 2); % Size: M x 1

% Find the worst outlier index
[max_score, bad_ensemble_idx] = max(anomaly_scores);
fprintf('The most abnormal trajectory is Index: %d with a score of %f\n', bad_ensemble_idx, max_score);


%% Instead of plotting the raw drift, color the trajectories by their total Anomaly Score!
spatial_anomaly_vector = zeros(MT_snap, 1);
for i = 1:M
    row_range = ((i-1)*steps_per_traj + 1) : (i*steps_per_traj);
    spatial_anomaly_vector(row_range) = anomaly_scores(i);
end

% Plot it! The "abnormal" trajectory will vividly light up in bright yellow/red
figure;
scatter3(Ya_all(1,:), Ya_all(2,:), Ya_all(3,:), 15, spatial_anomaly_vector, 'filled');
title('Lorenz Attractor Colored by Trajectory Anomaly Score');
colormap(hot); colorbar; view(45, 20); grid on;


% --- Bar Plot for Ensemble Index Inspection ---
figure('Color', 'w', 'Position', [200, 200, 600, 400]);

% 1. Plot all bars as standard blue first
hBar = bar(1:M, anomaly_scores, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'none'); 
hold on;

% 2. Highlight ONLY the abnormal trajectory in bright red
[max_score, bad_ensemble_idx] = max(anomaly_scores);
bar(bad_ensemble_idx, max_score, 'FaceColor', [0.9, 0.2, 0.2], 'EdgeColor', 'none');

% 3. Add a baseline or threshold text label over the bad bar
text(bad_ensemble_idx, max_score * 1.05, ['\leftarrow Outlier Trajectory #' num2str(bad_ensemble_idx)], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'Color', [0.9, 0.2, 0.2]);

% Formatting
xlabel('Ensemble Trajectory Index (i)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Anomaly Score (L_2 Drift Metric)', 'FontSize', 12, 'FontWeight', 'bold');
title('Dual Koopman Ensemble Anomaly Detection', 'FontSize', 13, 'FontWeight', 'bold');
set(gca, 'XTick', 1:M); % Ensure every trajectory index is labeled on X-axis
grid on;

%% --- 1. Compute Drift Vector Baseline ---
Ya_mean = zeros(3, steps_per_traj);
for m = 1:M
    Ya_mean = Ya_mean + Ya_ens{m};
end
Ya_mean = Ya_mean / M;

Ya_true_ens = cell(M, 1);
for m = 1:M
    Ya_true_ens{m} = Ya_ens{m} - Ya_mean;
end
Ya_true = cell2mat(Ya_true_ens.');

%% --- 2. Center the KEFs to isolate Drift ---
mean_KEFs = mean(KEFs, 1);
centered_KEFs_all = KEFs - mean_KEFs;

%% --- 3. Sort Modes by Drift Energy ---
total_modes_available = size(KEFs, 2);
mode_drift_power = zeros(total_modes_available, 1);
for j = 1:total_modes_available
    Y_mode_j = KMs(:, j) * Lambda_mat(j, j) * centered_KEFs_all(:, j).';
    mode_drift_power(j) = norm(Y_mode_j, 'fro');
end
[~, sort_drift_idx] = sort(mode_drift_power, 'descend');

%% --- 4. Compute Errors Per Trajectory & Normalize to 100% ---
max_J = min(90, total_modes_available); 
trajectory_errors_raw = zeros(M, max_J);

for j_test = 1:max_J
    active_modes = sort_drift_idx(1:j_test);
    Y_surrogate = KMs(:, active_modes) * Lambda_mat(active_modes, active_modes) * KEFs(:, active_modes).';
    Y_surr_tensor = reshape(Y_surrogate, 3, steps_per_traj, M);
    Y_surr_mean = mean(Y_surr_tensor, 3);
    
    for m = 1:M
        true_drift_m = Ya_true_ens{m};
        surr_drift_m = Y_surr_tensor(:,:,m) - Y_surr_mean;
        trajectory_errors_raw(m, j_test) = norm(true_drift_m - surr_drift_m, 'fro') / norm(true_drift_m, 'fro');
    end
end

% --- CRITICAL: Row-by-Row Normalization to force J=1 to be exactly 100% ---
trajectory_errors = zeros(M, max_J);
for m = 1:M
    % Divide each trajectory's profile by its own error at J=1, then scale to %
    trajectory_errors(m, :) = (trajectory_errors_raw(m, :) / trajectory_errors_raw(m, 1)) * 100;
end

% Pick standard clean columns to display for both plots
j_select = 1:5:max_J; 

%% --- 5. Plotting the Normalized Box Plot Evolution ---
figure('Color', 'w', 'Position', [100, 100, 700, 400]);
boxplot(trajectory_errors(:, j_select), j_select, 'Colors', [0.1 0.5 0.8]);
xlabel('Number of Drift-Sorted Koopman Modes (J)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Relative Drift Error Gap (%)', 'FontSize', 11, 'FontWeight', 'bold');
title('Box Plot: Ensemble Error Convergence Matrix (Normalized to 100%)', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0 105]);
grid on;

%% --- 6. Plotting the Native Violin Plot ---
figure('Color', 'w', 'Position', [150, 150, 750, 450]);
hold on;

num_vils = length(j_select);
x_spacing = 1:num_vils; % Clean spacing for plot positioning
width_scale = 0.4;     % Width of violin bodies

for i = 1:num_vils
    j_val = j_select(i);
    data_col = trajectory_errors(:, j_val);
    
    % Use Kernel Density Estimation to draw the distribution shape
    [f, y_grid] = ksdensity(data_col, 'NumPoints', 200);
    
    % Normalize density width to keep it tidy
    f = f / max(f) * width_scale;
    
    % Draw left and right patches to make it symmetric
    fill([x_spacing(i)-f, fliplr(x_spacing(i)+f)], [y_grid, fliplr(y_grid)], ...
         [0.2 0.6 0.8], 'FaceAlpha', 0.4, 'EdgeColor', [0.1 0.4 0.6], 'LineWidth', 1.5);
     
    % Overlay Mean and Median Indicators inside the violin
    plot(x_spacing(i), median(data_col), 'r*', 'MarkerSize', 8);
    plot([x_spacing(i)-width_scale/2, x_spacing(i)+width_scale/2], [mean(data_col), mean(data_col)], 'k-', 'LineWidth', 1.5);
end

% Set custom axis options to show true J numbers instead of loop indexes
set(gca, 'XTick', x_spacing, 'XTickLabel', j_select);
xlabel('Number of Drift-Sorted Koopman Modes (J)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Relative Drift Error Gap (%)', 'FontSize', 12, 'FontWeight', 'bold');
title('Violin Plot: Density and Spread of Ensemble Drift Errors', 'FontSize', 13, 'FontWeight', 'bold');
ylim([0 105]);
legend('Error Distribution Density', 'Median Indicator', 'Mean Indicator', 'Location', 'best');
grid on;

%%  --- 3. CRITICAL STEP: Sort Modes by Drift Energy ---
% We find which modes actually capture the deviation/drift matrix
% Compute a "Drift Power" score for each available Koopman mode
total_modes_available = size(KEFs, 2);
mode_drift_power = zeros(total_modes_available, 1);

for j = 1:total_modes_available
    % Single mode drift reconstruction
    Y_mode_j = KMs(:, j) * Lambda_mat(j, j) * centered_KEFs_all(:, j).';
    % Magnitude of drift captured by this specific mode
    mode_drift_power(j) = norm(Y_mode_j, 'fro');
end

% Sort indices based on maximum drift power (highest contribution first)
[~, sort_drift_idx] = sort(mode_drift_power, 'descend');

% --- 4. The Correct Truncation Loop (Using Sorted Indices) ---
max_J = min(30, total_modes_available); 
error_gap_vs_J = zeros(max_J, 1);

for j_test = 1:max_J
    % Pull the top j_test sorted modes that care about DRIFT
    active_modes = sort_drift_idx(1:j_test);
    
    % Reconstruct using only these drift-dominant modes
    Y_surrogate_drift = KMs(:, active_modes) * ...
                         Lambda_mat(active_modes, active_modes) * ...
                         centered_KEFs_all(:, active_modes).';
    
    % Compare drift vs drift
    error_gap_vs_J(j_test) = norm(Ya_true - Y_surrogate_drift, 'fro') / norm(Ya_true, 'fro');
end

% --- 5. Plotting the Decreasing Error Curve ---
figure('Color', 'w', 'Position', [100, 100, 550, 400]);
plot(1:max_J, error_gap_vs_J * 100, '-o', 'LineWidth', 2, 'Color', [0.1 0.5 0.8], 'MarkerFaceColor', [0.1 0.5 0.8]);
grid on;
xlabel('Number of Drift-Sorted Koopman Modes (J)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Surrogate Truncation Error Gap (%)', 'FontSize', 11, 'FontWeight', 'bold');
title('Surrogate Convergence via Drift-Targeted Ordering', 'FontSize', 12, 'FontWeight', 'bold');


% --- Compute the 2D Error Grid: Ensemble Index vs. Subspace Dimension J ---
max_J = min(30, total_modes_available);
ensemble_J_error_grid = zeros(M, max_J); % Rows = Ensembles, Cols = J

for j_test = 1:max_J
    active_modes = sort_drift_idx(1:j_test);
    
    % Reconstruct surrogate using the top j_test drift-sorted modes
    Y_surrogate = KMs(:, active_modes) * Lambda_mat(active_modes, active_modes) * KEFs(:, active_modes).';
    Y_surr_tensor = reshape(Y_surrogate, 3, steps_per_traj, M);
    Y_surr_mean = mean(Y_surr_tensor, 3);
    
    % Store error for each individual ensemble index
    for m = 1:M
        true_drift_m = Ya_true_ens{m};
        surr_drift_m = Y_surr_tensor(:,:,m) - Y_surr_mean;
        ensemble_J_error_grid(m, j_test) = norm(true_drift_m - surr_drift_m, 'fro') / norm(true_drift_m, 'fro') * 100;
    end
end

% --- Plotting the Ensemble Index Heatmap ---
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
imagesc(1:max_J, 1:M, ensemble_J_error_grid);
colormap(jet); % High error shows up as hot red, low error as cool blue
colorbar;
caxis([0 100]); % Set percentage limits

xlabel('Subspace Dimension (Number of Sorted Modes J)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Ensemble Index (m = 1 to M)', 'FontSize', 12, 'FontWeight', 'bold');
title('Universal Low-Order Convergence Matrix across Ensemble Indices', 'FontSize', 13, 'FontWeight', 'bold');
set(gca, 'YTick', 1:5:M); % Label ensemble indices every 5 ticks
grid on;

%% Trade-oof 

% Simulated Data for Demonstration (Replace with your actual loop results)
N_dict_range = [10, 20, 30, 50, 75, 100, 150, 200];
N_steps = T_snap;  % N time steps
M_ensembles = 10; % M trajectories

% 1. Compute metrics across N_dict
% Reconstruction loss (drops as N_dict increases)
loss_values = [0.45, 0.22, 0.12, 0.05, 0.03, 0.015, 0.009, 0.007]; 
% Computational Complexity Score proportional to (N * M * N_dict)
cost_values = (N_steps * M_ensembles) .* N_dict_range * 1e-3; 

% Choose your selected "optimal" point to highlight
selected_N_dict = 50;
idx_select = find(N_dict_range == selected_N_dict);

% --- START PLOT ---
figure('Color', 'w', 'Position', [100, 100, 650, 450]);
hold on;

% Left Axis: Loss/Error
yyaxis left
plot(N_dict_range, loss_values, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
ylabel('Reconstruction Loss / Error', 'FontSize', 11);
ax = gca; ax.YColor = 'b';

% Right Axis: Cost (N * M * N_dict footprint)
yyaxis right
plot(N_dict_range, cost_values, '-s', 'LineWidth', 2, 'MarkerFaceColor', 'r');
ylabel('Computational Footprint (N \times M \times N_{dict})', 'FontSize', 11);
ax = gca; ax.YColor = 'r';

% 2. Add the "Trade-off Cross-Line" at the selected point
% Vertical line marking the chosen N_dict
xline(selected_N_dict, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);

% Horizontal line for Left Axis (Loss at selected point)
yyaxis left
yline(loss_values(idx_select), ':', 'Color', 'b', 'LineWidth', 1.2);
plot(selected_N_dict, loss_values(idx_select), 'ok', 'MarkerSize', 10, 'LineWidth', 2); % highlight dot

% Horizontal line for Right Axis (Cost at selected point)
yyaxis right
yline(cost_values(idx_select), ':', 'Color', 'r', 'LineWidth', 1.2);
plot(selected_N_dict, cost_values(idx_select), 'sk', 'MarkerSize', 10, 'LineWidth', 2); % highlight square

% Formatting
xlabel('Number of Dictionary Functions (N_{dict})', 'FontSize', 12, 'FontWeight', 'bold');
title('Dual Koopman Trade-off & Capacity Selection', 'FontSize', 13, 'FontWeight', 'bold');
grid on;

%% -----------
figure('Color', 'w', 'Position', [150, 150, 600, 450]);
% Plot the frontier curve
plot(cost_values, loss_values, '-k', 'LineWidth', 1.5); hold on;

% Scatter points colored by their N_dict size to track capacity
scatter(cost_values, loss_values, 80, N_dict_range, 'filled', 'MarkerEdgeColor', 'k');
colormap(parula); colorbar;
hColor = colorbar; title(hColor, 'N_{dict}');

% Draw cross lines intersecting at your selected optimal operational point
xline(cost_values(idx_select), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
yline(loss_values(idx_select), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

% Annotate the selected point
text(cost_values(idx_select)*1.1, loss_values(idx_select)*1.3, ...
    ['\leftarrow Selected Optimal (N_{dict}=' num2str(selected_N_dict) ')'], ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');

xlabel('Computational Cost Factor (N \times M \times N_{dict})', 'FontSize', 11);
ylabel('Ensemble Reconstruction Loss', 'FontSize', 11);
title('Pareto Frontier: Accuracy vs. Computation Complexity', 'FontSize', 12, 'FontWeight', 'bold');
grid on;


%% 7. VISUALIZATION 3: Eigenfunction Spatial Slices
%%Build a dense 2D domain evaluation map at Z = 20
[gridX, gridY] = meshgrid(linspace(-20, 20, 100), linspace(-20, 20, 100));
X_test = [gridX(:)'; gridY(:)'; 20 * ones(1, 100*100)];

%%Calculate distance correlations using the native script's kernel function
G_test = kernel_f(X_test, Xa_all);

% Map the continuous slice to the first non-trivial eigenfunction track
% G_test is [10000 x 2400], G1 is [2400 x 2400] (or use pseudo-inverse logic)
% The mathematically clean projection for grid eigenfunctions:
Psi_grid = G_test * (G1 \ PX_trunc) * V_coeff(:, 1);

figure;
contourf(gridX, gridY, reshape(real(Psi_grid), 100, 100), 25, 'LineColor', 'none');
colormap(jet); colorbar;
title('Koopman Eigenfunction Spatial Map Slice (Z = 20)');
xlabel('X'); ylabel('Y');