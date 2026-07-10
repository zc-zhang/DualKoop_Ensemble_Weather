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
M=10;  % ensmeble
T_snap=size(Xa_all,2)/M;

idxM=(1:M);

N_dict = 450; % Choose your dictionary size (e.g., 150 features)
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
    %[V_coeff, Lambda, W2] = eig(K_star, G);   % G == eye(N_dict), so this ≈ eig(K_star)
    
     [V_coeff, Lambda] = eig(K_star, G);   % G == eye(N_dict), so this ≈ eig(K_star)
    KEs = diag(Lambda);
    KEFs = PX * V_coeff;


    %% Option 1: we nned W2m while it migbt be unstbale or ill-condition
   %RES1 = abs(sqrt(real(diag(W2'*L*W2)./diag(W2'*W2)-abs(KEs).^2)));

 %% Option 2: this exactly what Coolbrook did
  denominators = sum(abs(V_coeff).^2, 1).'; 
  numerators = real(sum(conj(V_coeff) .* (L * V_coeff), 1)).';
  RES = sqrt(max(0, (numerators ./ denominators) - abs(KEs).^2));


   % KMs     = Xa_all * pinv(KEFs).';
   KMs = Xa_all * pinv(KEFs.'); 
   % Ya_pred = KMs * Lambda * KEFs.';
  Ya_pred = KMs * diag(KEs) * KEFs.';


%%Calculate and display tracking error
reconstruction_error = norm(Ya_all - real(Ya_pred), 'fro') / norm(Ya_all, 'fro');
fprintf('=======================================\n');
fprintf('Trajectory Reconstruction Error: %.4f%%\n', reconstruction_error * 100);
fprintf('=======================================\n');



%%  VISUALIZATION 1:
figure('Position', [100, 100, 950, 420]);

%%Left Panel: True Dynamics
subplot(1, 2, 1);
scatter3(Ya_all(1,:), Ya_all(2,:), Ya_all(3,:), 12, Ya_all(3,:), 'filled');
title('First three elements (Ground Truth)');
xlabel('X'); ylabel('Y'); zlabel('Z');
view(45, 20); grid on; colorbar;

%%Right Panel: Reconstructed Dynamics
subplot(1, 2, 2);
scatter3(Ya_pred(1,:), Ya_pred(2,:), Ya_pred(3,:), 12, Ya_pred(3,:), 'filled');
%scatter3(Ya_all(1,:), Ya_all(2,:), Ya_all(3,:), 12, real(Ya_pred(3,:)), 'filled');
title(['Reconstructed three elements(Dual Koopman, N_d = ' num2str(N_dict) ')']);
xlabel('X'); ylabel('Y'); zlabel('Z');
view(45, 20); grid on; colorbar;

%% VISUALIZATION : Discrete Koopman Eigenvalues Spectrum
figure;
scatter(real(KEs),imag(KEs),300,RES,'.','LineWidth',1);
hold on
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
axis equal
axis([-1.15,1.15,-1.15,1.15])
clim([0,1])
load('cmap.mat')
colormap(jet); colorbar
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
%title(sprintf('Residuals ($M=%d$)',M),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;


%%   EIGVALS
figure;
theta = linspace(0, 2*pi, 100);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1.5); hold on;
scatter(real(KEs), imag(KEs), 45, 'r', 'filled');
title('Koopman Eigenvalue Distribution');
xlabel('Real(\lambda)'); ylabel('Imag(\lambda)');
grid on; axis equal;
legend('Unit Circle', 'Calculated KEs');


figure
loglog([0.001,1],[0.001,1],'k','linewidth',2)
hold on
loglog(sqrt(abs(abs(Lambda).^2-1)),RES,'b.','markersize',20)
xlabel('$\sqrt{|1-|\lambda|^2|}$','interpreter','latex','fontsize',18)
ylabel('residual','interpreter','latex','fontsize',18)
title(sprintf('Residuals ($N= %d$)',M),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;






%% some specific index of some location
indx_spec =1521;
figure; plot(Ya_all(1521,:),'b.-','LineWidth',1.2);
hold on; 
plot(real(Ya_pred(1521,:)),'r-','LineWidth',1.1); 




indx_spec = 1521;
figure; 

% 1. Create the shaded background regions FIRST so the lines sit on top
hold on;
colors = num2cell(lines(M), 2); % Generate M distinct colors for shading

for i = 1:M
    % Calculate the start and end point for this ensemble
    start_pt = (i-1)*T_snap + 1;
    end_pt = i*T_snap;
    
    % Draw a shaded background block for the current ensemble
    % (Using alpha to make it a soft, transparent pastel background)
    xr = xregion(start_pt, end_pt, 'FaceColor', colors{i}, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    
    % Add a text label near the top of each shaded block
    text(start_pt + T_snap/2, max(real(Ya_all(indx_spec,:))) * 0.9, ...
        sprintf('Ens# %d', i), 'HorizontalAlignment', 'center', ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', colors{i}*0.7);
end

% 2. Plot your actual data lines
p1 = plot(Ya_all(indx_spec,:), 'b.-', 'LineWidth', 1.2, 'DisplayName', 'True');
p2 = plot(real(Ya_pred(indx_spec,:)), 'r-', 'LineWidth', 1.1, 'DisplayName', 'Dual Koopman Pred');

% 3. Clean up the plot aesthetics
grid on;
xlim([1, size(Ya_all, 2)]);
xlabel('Total Snapshot Index (All Ensembles)', 'FontSize', 12);
ylabel(sprintf('State Value at Index %d', indx_spec), 'FontSize', 12);
title(sprintf('Ensemble Trajectory Comparison (State Variable %d)', indx_spec), 'FontSize', 14);
legend([p1, p2], 'Location', 'best');
ax = gca; ax.FontSize = 12;








%% Spatial modes
%% --- 1. User Settings and Grid Configurations ---
% FIX: Rows correspond to Z (97) and Columns correspond to Y (40)
grid_size = [97, 40];  % [Rows (z), Columns (y)] -> 97 * 40 = 3880

% Specify your target coordinates directly based on the true vectors
highlight_row = 25;    % Position along the vertical z-axis (1 to 97)
highlight_col = 26;    % Position along the horizontal y-axis (1 to 40)

% Let MATLAB natively calculate the exact 1D index using the corrected grid size
highlight_index = sub2ind(grid_size, highlight_row, highlight_col);

% Choose your sorting criteria: 'RES' or 'Magnitude'
sorting_criteria = 'Magnitude'; 

%% --- 2. Automated Dual Koopman Sorting Mechanism ---
if strcmp(sorting_criteria, 'RES')
    [~, sort_idx] = sort(RES, 'ascend');
    plot_title_prefix = 'True Mode (Low RES) ';
else
    %[~, sort_idx] = sort(abs(abs(KEs) - 1), 'ascend'); 
     [~, sort_idx] = sort(abs(KEs), 'ascend'); 
    plot_title_prefix = 'Persistent Mode ';
end

KMs_sorted = KMs(:, sort_idx);

%% --- 3. Figure Layout and Spatial Visualization ---
figure;
numRows = 2;
numCols = 3;
maxTiles = numRows * numCols;
ttt = tiledlayout(numRows, numCols, 'Padding', 'compact', 'TileSpacing', 'compact');

modesToPlot = 1:maxTiles; 

for k = 1:length(modesToPlot)
    mode_index_in_sorted = modesToPlot(k);
    
    % Extract the spatial mode column vector (Size: 3880 x 1)
    KMode_i = abs(KMs_sorted(:, mode_index_in_sorted)); 
    
    % Reshape perfectly back to your [97, 40] physical grid layout
    mode_i_reshaped = reshape(KMode_i, grid_size); 
    
    ax = nexttile;
    
    % FIX: imagesc expects (X_vector, Y_vector, Matrix) 
    % where Matrix dimensions must be length(Y_vector) x length(X_vector).
    % Since data.y is 40x1 and data.z is 97x1, and mode_i_reshaped is 97x40:
    imagesc(data.y, data.z, mode_i_reshaped); 
    
    colormap(ax, brighten(redblueTecplot(21), -0.55));
    colorbar;
    
    % Maintain physical coordinates tracking orientation (z going upwards)
    axis xy;
    xlim([0 2e4]);
    ylim([0 2e4]);
    
    set(gca, 'FontSize', 17.6);  
    xlabel("y (m)", 'FontSize', 17.6); % y is horizontal
    ylabel("z (m)", 'FontSize', 17.6); % z is vertical
    
    orig_idx = sort_idx(mode_index_in_sorted);
    title(sprintf('%s#%d (\\lambda=%.2f)', plot_title_prefix, k, abs(KEs(orig_idx))), 'FontSize', 14);

    %% --- 4. Direct Coordinate Marker Overlay ---
    % Grab coordinate positions straight from your corrected indices
    highlight_y = data.y(highlight_col); 
    highlight_z = data.z(highlight_row); 
    
    hold on;
    plot(highlight_y, highlight_z, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow');
    hold off;
end

set(gcf, 'Renderer', 'painters');






%% --- Simultaneous Spatial Comparison Across Ensembles ---
grid_size = [97, 40];  % [Rows (z), Columns (y)]

% Target settings
target_time = 25;       % Choose time step: 5 (initial stage) or 45 (late stage)
ensembles_to_check = 1:6; % Simultaneously look at the first 6 ensemble members

figure;
t_layout = tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t_layout, sprintf('Spatial Field Comparison at Time %d (N_{dict} = %d)', ...
    target_time, size(KMs, 2)), 'FontSize', 16, 'FontWeight', 'bold');

for k = 1:length(ensembles_to_check)
    current_ens = ensembles_to_check(k);
    
    % 1. Calculate the exact snapshot column index for this specific ensemble
    snapshot_idx = (current_ens - 1) * T_snap + target_time;
    
    % 2. Reconstruct the spatial field vector for this snapshot
    snapshot_eigenfunctions = KEFs(snapshot_idx, :); 
    spatial_field_vector = real(KMs * diag(KEs) * snapshot_eigenfunctions.'); 
    
    % 3. Reshape back to the true physical [97 x 40] grid
    spatial_field_reshaped = reshape(spatial_field_vector, grid_size);
    
    % Activate the next tile in our layout
    ax = nexttile;
    imagesc(data.y, data.z, spatial_field_reshaped);
    colormap(ax, brighten(redblueTecplot(21), -0.55));
    colorbar;
    axis xy;
    xlim([0 2e4]); ylim([0 2e4]);
    
    set(gca, 'FontSize', 12);
    xlabel("y (m)", 'FontSize', 11);
    ylabel("z (m)", 'FontSize', 11);
    title(sprintf('Ensemble Member %d', current_ens), 'FontSize', 13);
    
    %% --- Coordinate Marker Overlay ---
    highlight_row = 25; highlight_col = 26;
    highlight_y = data.y(highlight_col); 
    highlight_z = data.z(highlight_row); 
    hold on;
    plot(highlight_y, highlight_z, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow');
    hold off;
end

set(gcf, 'Renderer', 'painters');









%% --- 2. Log-Scale Visualization: Bar Chart & Distribution Contrast ---
figure('Position', [100, 100, 950, 450]);

% --- LEFT: Logarithmic Bar Chart ---
subplot(1,2,1);
b = bar(1:M, ensemble_anomaly_bar, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none');
grid on; hold on;

% Highlight the abnormal index in bright red
bar(abnormal_member_idx, max_val, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none');

% CHANGE TO LOG SCALE
set(gca, 'YScale', 'log'); 

% Clean up the labels and layout
set(gca, 'XTick', 1:M, 'FontSize', 12);
xlabel('Ensemble Index (M = 1:10)', 'FontWeight', 'bold', 'FontSize', 13);
ylabel('Mean Raw Data Deviation Norm (Log Scale)', 'FontWeight', 'bold', 'FontSize', 13);
title('Log-Scale Outlier Index Identification', 'FontSize', 13);
legend('Normal Cluster', sprintf('Abnormal Member (Index %d)', abnormal_member_idx), 'Location', 'best');

% Adjust Y-limits slightly so the bars don't clip at the bottom of the log axis
min_val = min(ensemble_anomaly_bar);
ylim([min_val * 0.5, max_val * 2]);

% --- RIGHT: Logarithmic Box & Swarm Plot ---
subplot(1,2,2);
normal_group_data = ensemble_anomaly_bar(1:M ~= abnormal_member_idx);

% Plot the boxplot on a log scale
boxplot([normal_group_data; max_val], [ones(9,1); 2], 'Labels', {'Normal Group (9)', 'Outlier (1)'});
grid on; hold on;

% Superimpose the raw data points as a distribution overlay
swarmchart(ones(9,1), normal_group_data, 50, [0.3 0.5 0.8], 'filled', 'XJitter', 'density');
swarmchart(2, max_val, 90, [0.8 0.2 0.2], 'filled');

set(gca, 'YScale', 'log');
set(gca, 'FontSize', 12);
ylabel('Mean Raw Data Deviation Norm (Log Scale)', 'FontWeight', 'bold', 'FontSize', 13);
title('Log-Scale Statistical Separation', 'FontSize', 13);
ylim([min_val * 0.5, max_val * 2]);


%%   % Additional RESIDUAL TESTINg for continous spetra
% x_pts=-1.2:0.02:1.2;    y_pts=-0.02:0.02:1.2;
% z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);		% complex points where we compute pseudospectra
% RES0 = KoopPseudoSpecQR(PX,PY,1/N_dict,z_pts);
% RES0=reshape(RES0,length(y_pts),length(x_pts));
% 
% RES1 = KoopPseudoSpec(double(G),double(K_star),double(L),z_pts,'Parallel','off');	% compute pseudospectra
% RES1=reshape(RES1,length(y_pts),length(x_pts));
% 
% %% Plot pseudospectra
% figure
% hold on
% v=(10.^(-10:0.2:0));
% contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES0)),log10(v));
% hold on
% contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES0)),log10(v));
% cbh=colorbar;
% cbh.Ticks=log10(10.^(-4:1:0));
% cbh.TickLabels=10.^(-4:1:0);
% clim([-4,0]);
% reset(gcf)
% set(gca,'YDir','normal')
% colormap gray
% axis equal;
% 
% title(sprintf('Naive Residual ($N_dict=%d$)',M),'interpreter','latex','fontsize',18)
% xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
% ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
% 
% ax=gca; ax.FontSize=18; axis equal tight;   axis([x_pts(1),x_pts(end),-y_pts(end),y_pts(end)])
% hold on
% plot(real(KEs),imag(KEs),'.r','markersize',12);
% box on
% 
% 
% 
% 
% 
% figure
% hold on
% v=(10.^(-10:0.2:0));
% contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES)),log10(v));
% hold on
% contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES)),log10(v));
% cbh=colorbar;
% cbh.Ticks=log10(10.^(-4:1:0));
% cbh.TickLabels=10.^(-4:1:0);
% clim([-4,0]);
% reset(gcf)
% set(gca,'YDir','normal')
% colormap gray
% axis equal;
% 
% title(sprintf('Pseudospectrum ($N_d=%d$)',N_dict),'interpreter','latex','fontsize',18)
% xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
% ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
% 
% ax=gca; ax.FontSize=18; axis equal tight;   axis([x_pts(1),x_pts(end),-y_pts(end),y_pts(end)])
% hold on
% plot(real(KEs),imag(KEs),'.r','markersize',12);
% box on


























%% Additional RESIDUAL compute Actua;lyy we have doen in R
%% --- 1. Compute the Residuals (Your Current Block) ---
% denominators = sum(abs(V_coeff).^2, 1).'; 
% numerators = real(sum(conj(V_coeff) .* (L * V_coeff), 1)).';
% RES = sqrt(max(0, (numerators ./ denominators) - abs(KEs).^2));

% RES = zeros(N_dict, 1);
% for j = 1:N_dict
%     v   = V_coeff(:, j);
%     lam = KEs(j);
%     % Full residual: ||Ug - λg||² / ||g||²
%     res_num = real(v' * (L - conj(lam)*K_star - lam*K_star' + abs(lam)^2*G) * v);
%     res_den = real(v' * G * v);
%     RES(j)  = sqrt(max(0, res_num / res_den));
% end
% %[~, sort_res_idx] = sort(RES, 'ascend');  % small residual = genuine eigenvalue
% 
% %% --- 2. Sort by Residual in ASCENDING Order (Best Physics First) ---
% [sorted_RES, sort_res_idx] = sort(RES, 'ascend');
% 
% %% --- 3. Evaluate Reconstruction Across Your Truncation Set ---
% %%  Compute Trajectory Errors ---
% N_dict_trunc_set = [10, 20, 40, 60, 90, 300]; 
% num_tests = length(N_dict_trunc_set);
% Ya_all_tensor = reshape(real(Ya_all), p, T_snap, M);
% error_data = zeros(M, num_tests);
% 
% % Initialize cell array to store the 3D tensor for each test size
% Ya_pred_N_ditcs = cell(num_tests, 1);
% 
% for i = 1:num_tests
%     active_indices = sort_res_idx(1:N_dict_trunc_set(i));
% 
%     % Quick Reconstruction Slicing
%     V_trunc = V_coeff(:, active_indices);
%     L_trunc = diag(KEs(active_indices));
%     Ya_pred_t = (Xa_all * pinv(PX * V_trunc).') * L_trunc * (PX * V_trunc).';
%     Ya_pred_tensor = reshape(real(Ya_pred_t), 3, T_snap, M);
% 
%     Ya_pred_N_ditcs{i} = Ya_pred_tensor;
%     % Track individual trajectory errors
%     for m = 1:M
%         error_data(m, i) = (norm(Ya_all_tensor(:,:,m) - Ya_pred_tensor(:,:,m), 'fro') / ...
%                             norm(Ya_all_tensor(:,:,m), 'fro')) * 100;
%     end
% end
% 
% %% --- 2. Plot Both Side-by-Side (Ultra Simple) ---
% figure('Color', 'w', 'Position', [100, 100, 1100, 500]);
% x_labels = arrayfun(@(x) ['N=', num2str(x)], N_dict_trunc_set, 'UniformOutput', false);
% 
% % LEFT PANEL: Simple Standard Boxplot
% %subplot(1, 2, 1);
% boxplot(error_data, 'Labels', x_labels, 'Widths', 0.5);
% title('Boxplot Error Distribution', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Relative Reconstruction Error (%)', 'FontSize', 11);
% xlabel('Active RKHS Features (N)', 'FontSize', 11);
% ylim([0, 105]);
% grid on;
% 
% 
% 
% %% --- Define Which Loop Indices You Want to Plot ---
% idx_Ndict1 = 5; % Pulls the 2nd dictionary size from your set (e.g., N = 20)
% idx_Ndict2 = 6; % Pulls the 6th dictionary size from your set (e.g., N = 100)
% 
% %% --- Widescreen Plotting via Direct Cell Extraction ---
% track_colors = jet(M); 
% 
% figure('Color', 'w', 'Position', [100, 100, 1200, 550]);
% 
% % =========================================================================
% % SUBPLOT 1: Ground Truth
% % =========================================================================
% subplot(1, 3, 1);
% hold on;
% for m = 1:M
%     % 1. Plot the trajectory line
%     plot3(Ya_all_tensor(1,:,m), Ya_all_tensor(2,:,m), Ya_all_tensor(3,:,m), ...
%           'Color', [track_colors(m, :), 0.4], 'LineWidth', 1.5); 
%     % 2. Highlight the Initial Condition (Time Step 1)
%     plot3(Ya_all_tensor(1,1,m), Ya_all_tensor(2,1,m), Ya_all_tensor(3,1,m), ...
%           'o', 'MarkerFaceColor', track_colors(m, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
% end
% 
% % --- ADDED: Compute and Overlay True Ensemble Mean ---
% Ya_mean_true = mean(Ya_all_tensor, 3); % Size: [3 x T_snap]
% plot3(Ya_mean_true(1,:), Ya_mean_true(2,:), Ya_mean_true(3,:), 'k--o', ...
%       'LineWidth', 2.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k');
% 
% title('True Ensemble (Ground Truth)', 'FontSize', 11, 'FontWeight', 'bold');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% view(45, 20); grid on;
% hold off;
% 
% % =========================================================================
% % SUBPLOT 2: First Selected Dictionary Truncation
% % =========================================================================
% subplot(1, 3, 2);
% hold on;
% for m = 1:M
%     % 1. Plot the reconstructed trajectory line
%     plot3(Ya_pred_N_ditcs{idx_Ndict1}(1,:,m), Ya_pred_N_ditcs{idx_Ndict1}(2,:,m), Ya_pred_N_ditcs{idx_Ndict1}(3,:,m), ...
%           'Color', [track_colors(m, :), 0.4], 'LineWidth', 1.5); 
%     % 2. Highlight the Initial Condition (Time Step 1)
%     plot3(Ya_pred_N_ditcs{idx_Ndict1}(1,1,m), Ya_pred_N_ditcs{idx_Ndict1}(2,1,m), Ya_pred_N_ditcs{idx_Ndict1}(3,1,m), ...
%           'o', 'MarkerFaceColor', track_colors(m, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
% end
% 
% % --- ADDED: Compute and Overlay Reconstructed Mean (Truncated) ---
% Ya_mean_pred1 = mean(Ya_pred_N_ditcs{idx_Ndict1}, 3);
% plot3(Ya_mean_pred1(1,:), Ya_mean_pred1(2,:), Ya_mean_pred1(3,:), 'k--o', ...
%       'LineWidth', 2.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k');
% % -----
% title(['Reconstructed (N = ' num2str(N_dict_trunc_set(idx_Ndict1)) ')'], 'FontSize', 11, 'FontWeight', 'bold');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% view(45, 20); grid on;
% hold off;
% 
% % =========================================================================
% % SUBPLOT 3: Second Selected Dictionary Truncation
% % =========================================================================
% subplot(1, 3, 3);
% hold on;
% for m = 1:M
%     % 1. Plot the reconstructed trajectory line
%     plot3(Ya_pred_N_ditcs{idx_Ndict2}(1,:,m), Ya_pred_N_ditcs{idx_Ndict2}(2,:,m), Ya_pred_N_ditcs{idx_Ndict2}(3,:,m), ...
%           'Color', [track_colors(m, :), 0.4], 'LineWidth', 1.5); 
%     % 2. Highlight the Initial Condition (Time Step 1)
%     plot3(Ya_pred_N_ditcs{idx_Ndict2}(1,1,m), Ya_pred_N_ditcs{idx_Ndict2}(2,1,m), Ya_pred_N_ditcs{idx_Ndict2}(3,1,m), ...
%           'o', 'MarkerFaceColor', track_colors(m, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
% end
% 
% % --- ADDED: Compute and Overlay Reconstructed Mean (Full Dictionary) ---
% Ya_mean_pred2 = mean(Ya_pred_N_ditcs{idx_Ndict2}, 3);
% plot3(Ya_mean_pred2(1,:), Ya_mean_pred2(2,:), Ya_mean_pred2(3,:), 'k--o', ...
%       'LineWidth', 2.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k');
% 
% title(['Reconstructed (N = ' num2str(N_dict_trunc_set(idx_Ndict2)) ')'], 'FontSize', 11, 'FontWeight', 'bold');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% view(45, 20); grid on;
% hold off;
