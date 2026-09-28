%% Ensemble data collection


%% 
%latN = 97;
%lonN = 40;

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


%% Ensemable MEan data handle 

% Parameters
%p = 3880;
%T_snap = 50;
%M = 10;

% Xa_EnMean = zeros(p, T_snap);
% Ya_EnMean = zeros(p, T_snap);
% 
% for i = 1:M
%     col_indices = (i-1)*T_snap + (1:T_snap);
%     Xa_EnMean = Xa_EnMean + Xa_all(:, col_indices);
%     Ya_EnMean = Ya_EnMean + Ya_all(:, col_indices);
% end
% EnMVfull_Xa = Xa_EnMean / M;
% EnMVfull_Ya = Ya_EnMean / M;
% 
% 
% % 
% % % 3. Save the resulting mean data
%  save(fullfile(baseDir, 'EnMVfull_Xa.mat'), 'EnMVfull_Xa', '-v7.3');
%  save(fullfile(baseDir, 'EnMVfull_Ya.mat'), 'EnMVfull_Ya', '-v7.3');
% 
% fprintf('Ensemble Mean matrices derived successfully.\n');
% fprintf('New dimensions: %d x %d\n', size(EnMVfull_Xa, 1), size(EnMVfull_Xa, 2));


%clc; clear all
%shapefile_path = 'D:\Susuki Lab\Testing_Code\data-weather\Data_250401\ne_10m_coastline\ne_10m_coastline.shp';
%S = shaperead(shapefile_path);



%% Load (EnVfull_all260708.mat)
Xa_all =EnVfull_Xa_36;
Ya_all =EnVfull_Ya_36;


%=================
% compute the ensemble mean 
Xa_raw_3D               = reshape(Ya_all, 3880, 36, 10);
Xa_mean         = mean(Xa_raw_3D, 3);  
Ya_raw_3D               = reshape(Ya_all, 3880, 36, 10);
Ya_mean         = mean(Ya_raw_3D, 3);  

%=======================
%Xa_EnMean =EnMVfull_Xa;
%Ya_EnMean=EnMVfull_Ya;

p=size(EnVfull_Ya_36,1);
M=10;  % ensmeble
T_snap=size(EnVfull_Ya_36,2)/M;

idxM=(1:M);

N_dict = 330; % Choose your dictionary size (e.g., 150 features)
fprintf('Running kernel_ResDMD with N = %d features...\n', N_dict);

%%Run the original code signature to get the clean truncated spaces
%[G, K_star, L, PX, PY, PSI_x, PSI_y, PSI_y2, G1, A1, kernel_f] = ...
   % kernel_ResDMD(Xa_all, Ya_all, 'type', 'Gaussian', 'N', N_dict);
[G, K_star, L, PX, PY, PSI_x0, PSI_y0, ~, G1, A1, UU, kernel_f] = kernel_ResDMD(...
    Xa_all, Ya_all, ...
    'type',    'Laplacian', ...
    'N',       N_dict, ...    %'cut_off', 1e-6, ...
    'Xb',      Xa_mean, ...  % Out-of-sample input (Ensemble Mean t=0)
    'Yb',      Ya_mean);     % Out-of-sample target (Ensemble Mean t=t_1)

% 3. COMPUTE THE DISCRETE FEATURE COMPONENTS HERE (FIXED)
fprintf('Computing KEs, KEFs, and KMs on the truncated space...\n');

%% ========================= NOn TRunction ============
    % No truncation needed — G, K_star, PX are already N_dict-sized
    %[V_coeff, Lambda, W2] = eig(K_star, G);   % G == eye(N_dict), so this ≈ eig(K_star)
    
     [V_coeff, Lambda, W_coeff] = eig(K_star);   % G == eye(N_dict), so this ≈ eig(K_star)
    KEs = diag(Lambda);
    KEs_dual = conj(KEs);           % vector — conjugate of KEs
    KEFs = PX * V_coeff;
    KEFs_dual = PX * W_coeff;        % φ_j-dual = φ_j(x_i) evaluated via dual — i.e. ⟨φ_j, κ_{x_i}⟩

    %% Option 1: we nned W2m while it migbt be unstbale or ill-condition
   %RES1 = abs(sqrt(real(diag(W2'*L*W2)./diag(W2'*W2)-abs(KEs).^2)));

 %% Option 2: this exactly what Coolbrook did
  denominators = sum(abs(V_coeff).^2, 1).'; 
  numerators = real(sum(conj(V_coeff) .* (L * V_coeff), 1)).';
  RES = sqrt(max(0, (numerators ./ denominators) - abs(KEs).^2));


   % KMs     = Xa_all * pinv(KEFs).';
   KMs = Xa_all * pinv(KEFs.');   %%  KMs = Xa_all / (KEFs.');
   % Ya_pred = KMs * Lambda * KEFs.';
  Ya_pred = KMs * diag(KEs) * KEFs.';

% Dual KMD 
% Eq. (2): y_t^i = Σ_j λ̄_j^t · coeffs(:,j) · KEFs_dual(i,j)
Ya_pred_dual = KMs * diag(KEs_dual) * KEFs_dual.';

%%Calculate and display tracking error
reconstruction_error = norm(Ya_all - real(Ya_pred), 'fro') / norm(Ya_all, 'fro');
fprintf('=======================================\n');
fprintf('Trajectory Reconstruction Error: %.4f%%\n', reconstruction_error * 100);
fprintf('=======================================\n');


%%Calculate and display tracking error
dual_recons_error = norm(Ya_all - real(Ya_pred_dual), 'fro') / norm(Ya_all, 'fro');
fprintf('=======================================\n');
fprintf('Dual Traj Reconstruction Error: %.4f%%\n', reconstruction_error * 100);
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
%scatter3(Ya_all(1,:), Ya_all(2,:), Ya_all(3,:), 12, real(Ya_pred(3,:)),
%'filled');
title(['Reconstructed three elements(Dual Koopman, N_d = ' num2str(N_dict) ')']);
xlabel('X'); ylabel('Y'); zlabel('Z');
view(45, 20); grid on; colorbar;

%% VISUALIZATION : Discrete Koopman Eigenvalues Spectrum
figure;
%scatter(real(KEs),imag(KEs),300,RES,'.','LineWidth',1);
scatter(real(KEs_dual),imag(KEs_dual),300,RES,'.','LineWidth',1);
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
ax=gca; ax.FontSize=18; box on;


%%   EIGVALS
figure;
theta = linspace(0, 2*pi, 100);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1.5); hold on;
%scatter(real(KEs), imag(KEs), 45, 'r', 'filled');
scatter(real(KEs_dual), imag(KEs_dual), 45, 'r', 'filled');
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
title(sprintf('Residuals (N= num2str(N_dict))',M),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;box on;






%% some specific index of some location
indx_spec = 2450; % 1521;indx_spec = 2450; % 1521;
figure; 
plot(Ya_all(indx_spec,:),'b.-','LineWidth',1.2);
hold on; 
plot(real(Ya_pred(indx_spec,:)),'r-','LineWidth',1.1); 
hold on;
plot(real(real(Ya_pred_dual(indx_spec,:))),'g-','LineWidth',1.1); 


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

%% Choose your sorting criteria: 'RES', 'Magnitude', or 'Periodic'
sorting_criteria = 'Magnitude'; 
dt = 1; % <--- Set your true Delta t simulation time step here

%% --- 2. Automated Dual Koopman Sorting Mechanism ---
if strcmp(sorting_criteria, 'RES')
    [~, sort_idx] = sort(RES, 'ascend');
    plot_title_prefix = 'True Mode (Low RES) ';
    
elseif strcmp(sorting_criteria, 'Magnitude')
    [~, sort_idx] = sort(abs(KEs), 'descend'); 
  %  [~, sort_idx] = sort(abs(KEs_dual), 'descend'); 
    plot_title_prefix = 'Dominant Mode ';
    
elseif strcmp(sorting_criteria, 'Periodic')
    % Direct mathematical evaluation of 2*pi / omega
    T_periodic = (2 * pi * dt) ./ abs(angle(KEs));
     %
     % T_periodic = (2 * pi * dt) ./ abs(angle(KEs_dual));
    % Sort directly by raw value (oscillations first, non-oscillating/Inf last)
    [~, sort_idx] = sort(T_periodic, 'descend'); 
    plot_title_prefix = 'Periodic Mode ';
end

% Perfectly synchronized re-ordering of all core arrays
KEs  = KEs(sort_idx);        
KMs  = KMs(:, sort_idx);     
KEFs = KEFs(:, sort_idx);

% KEs_dual  = KEs(sort_idx);        
% KMs  = KMs(:, sort_idx);     
% KEFs_dual = KEFs(:, sort_idx);
% 
%===
RES_sorted=RES(sort_idx);

%% --- 3. Figure Layout and Spatial Visualization ---
figure;
numRows = 3;
numCols = 3;
maxTiles = numRows * numCols;
ttt = tiledlayout(numRows, numCols, 'Padding', 'compact', 'TileSpacing', 'compact');

% Skip the complex conjugates to only plot unique spatial patterns!
modesToPlot = 1:2:17; 

for k = 1:length(modesToPlot)
    % This is the TRUE sorted index (1, 3, 5, 7, ..., 17)
    mode_index_in_sorted = modesToPlot(k);
    
    % Extract the spatial mode column vector (Size: 3880 x 1)
    KMode_i = abs(KMs(:, mode_index_in_sorted)); 
    
    % Reshape perfectly back to your physical grid layout
    mode_i_reshaped = reshape(KMode_i, grid_size); 
    
    % Automatically jumps to the next available slot in the 3x3 layout
    ax = nexttile;
    
    imagesc(data.y, data.z, mode_i_reshaped); 
    colormap(ax, brighten(redblueTecplot(21), -0.55));
    colorbar;
    
    % Maintain physical coordinates tracking orientation (z going upwards)
    axis xy;
    axis equal; 
    xlim([0 2e4]);
    ylim([0 2e4]);
    
    set(gca, 'FontSize', 15);  
    xlabel("y", 'FontSize', 10,'latex'); 
    ylabel("z", 'FontSize', 10,'latex'); 
    
   % 1. Extract the current complex eigenvalue
    lambda_j = KEs(mode_index_in_sorted);
   % lambda_j = KEs_dual(mode_index_in_sorted);
    real_part = real(lambda_j);
    imag_part = imag(lambda_j);
    
    % 2. Dynamic LaTeX String: %d automatically changes the subscript to 1, 3, 5...
    title_str = sprintf('$\\textrm{Mode }%d: \\lambda_{%d} = %.2f %+.2f i$', ...
        mode_index_in_sorted, mode_index_in_sorted, real_part, imag_part,'interpreter', 'latex', 'fontsize', 17);
    
    % 3. Plot with LaTeX Interpreter
    title(title_str, 'Interpreter', 'latex', 'FontSize', 10);
    
    %% --- 4. Direct Coordinate Marker Overlay ---
    highlight_y = data.y(highlight_col); 
    highlight_z = data.z(highlight_row); 
    % 
     hold on;
     plot(highlight_y, highlight_z, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow');
     hold off;
end
set(gcf, 'Renderer', 'painters');
    exportgraphics(gcf, sprintf('ensemble_mode_%d.pdf', j), ...
        'ContentType', 'vector', 'BackgroundColor', 'none')





%% Raw data on Accumulated Vorticity of Ensemble Mean
figure;
Ya_EnMean=Ya_mean;
% === Accumulation over time (sum instead of mean) ===
spatial_accumulated = reshape(sum(Ya_EnMean, 2), [97,40]);

% Automatically jumps to the next available slot in the 3x3 layout
ax = nexttile;
imagesc(data.y, data.z, spatial_accumulated);

colormap(ax, brighten(redblueTecplot(21), -0.55));
colorbar;

% Maintain physical coordinates
axis xy;
axis equal;
xlim([0 2e4]);
ylim([0 2e4]);
set(gca, 'FontSize', 11);
xlabel("y", 'FontSize', 11);
ylabel("z", 'FontSize', 11);
title(sprintf('Accumulated Vorticity of Ens. Mean'), 'FontSize', 10);

%% --- Coordinate Marker Overlay ---
highlight_row = 25; 
highlight_col = 26;
highlight_y = data.y(highlight_col);
highlight_z = data.z(highlight_row);

hold on;
plot(highlight_y, highlight_z, 'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'yellow');
hold off;

%% --- Simultaneous Spatial Comparison Across Ensembles ---
grid_size = [97, 40];  % [Rows (z), Columns (y)]

% Target settings
target_time = 25;       % Choose time step: 5 (initial stage) or 45 (late stage)
ensembles_to_check = 1:9; % Simultaneously look at the first 6 ensemble members

figure;

% Automatically calculate layout dimensions (e.g., 3 columns, and as many rows as needed)
num_plots = length(ensembles_to_check);
num_cols = 3; 
num_rows = ceil(num_plots / num_cols);

t_layout = tiledlayout(num_rows, num_cols, 'Padding', 'compact', 'TileSpacing', 'compact');
%t_layout = tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t_layout, sprintf('Spatial Field Comparison at Time %d (N_{dict} = %d)', ...
    target_time, size(KMs, 2)), 'FontSize', 16, 'FontWeight', 'bold');

for k = 1:length(ensembles_to_check)
    current_ens = ensembles_to_check(k);
    

    %% Option 1: use the KEF embeeding the time involved --\varphi_j(x_t^i)
    % % 1. Calculate the exact snapshot column index for this specific ensemble
     %snapshot_idx = (current_ens - 1) * T_snap + target_time;
    % 
    % % 2. Reconstruct the spatial field vector for this snapshot
    % snapshot_eigenfunctions = KEFs(snapshot_idx, :); % i.e., phi_t = KEFs(snapshot_idx, :); 
    %%%%% WRONG： spatial_field_vector = real(KMs * diag(KEs) * snapshot_eigenfunctions.'); 
    %spatial_field_vector = real(KMs * snapshot_eigenfunctions.'); 

    %% Option 2: use the time evolution by lambda_j^t<\varphi_j,\kappa_{x_0^i}>
    % 1. Find the index of the INITIAL condition (t = 1 or t = 0) for this ensemble
    % Assumes the first snapshot of each ensemble block in KEFs is the initial state
      initial_snapshot_idx = (current_ens - 1) * T_snap + 1; 
     phi_zero = KEFs(initial_snapshot_idx, :); % This is <\phi_j, \kappa_{x_0^i}>
    % 
    % % 2. Pure Koopman Evolution: Propagate forward using \lambda_j^(target_time)
    % % Note: If your system matches t=1 as initial, use (target_time - 1)
      t_power = target_time - 1; 
      phi_t = phi_zero .* (KEs(:).').^t_power; 
    % 
    % % 3. Reconstruct back into physical space via Koopman Modes (KMs)
    % % (Note: If KEs are already folded into KMs or phi_t handles the evolution, 
    % % you just multiply KMs by phi_t)
     spatial_field_vector = real(KMs * phi_t.');



    % 3. Reshape back to the true physical [97 x 40] grid
    spatial_field_reshaped = reshape(spatial_field_vector, grid_size);
    
    % Activate the next tile in our layout
    ax = nexttile;
    imagesc(data.y, data.z, spatial_field_reshaped);
    colormap(ax, brighten(redblueTecplot(21), -0.55));
    colorbar;
    axis xy;
    axis equal; % <-- Add this to stop the stretching!
    xlim([0 2e4]); ylim([0 2e4]);
    
    set(gca, 'FontSize', 12);
    xlabel("y", 'FontSize', 11);
    ylabel("z", 'FontSize', 11);
    title(sprintf('Ensemble # %d', current_ens), 'FontSize', 13);
    
    %% --- Coordinate Marker Overlay ---
    highlight_row = 25; highlight_col = 26;
    highlight_y = data.y(highlight_col); 
    highlight_z = data.z(highlight_row); 
    hold on;
    plot(highlight_y, highlight_z, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow');
    hold off;
end

set(gcf, 'Renderer', 'painters');


%% Option 1 emphaseiz
%% check y_t,p^i =<h_p,\kappa_{x_t^i}>
% Target settings
target_time = 25;         % Target snapshot time-step to visualize
ensembles_to_check = 1:9; % Look at the first 9 ensemble members simultaneously

grid_size = [97, 40];     % [Rows (z), Columns (y)] -> Total 3880 states
T_snap = 36;              % Your preset snapshot horizon length per ensemble

% --- CRITICAL SUBSPACE FIX ---
% To ensure the output operator works cleanly without cross-column scrambling:
r = min(150, size(KEFs, 2)); % Set an appropriate rank truncation r
KEs_r  = KEs(sort_idx(1:r));
KEFs_r = KEFs(:, sort_idx(1:r));
KMs_r  = EnVfull_Xa_36 / (KEFs_r.'); % Clean algebraic right-division for accurate spatial modes

figure;%('Name', 'Dual Koopman Spatial Diagnostics (Option 1)', 'Position', [100, 100, 1200, 850]);

% Automatically calculate layout dimensions for 3 columns
num_plots = length(ensembles_to_check);
num_cols = 3; 
num_rows = ceil(num_plots / num_cols);
t_layout = tiledlayout(num_rows, num_cols, 'Padding', 'compact', 'TileSpacing', 'compact');

% title(t_layout, sprintf('Dual Koopman y_t^i = C_h(\\kappa_{x_t^i}) at Time %d (Rank r = %d)', ...
%     target_time, r), 'FontSize', 16, 'FontWeight', 'bold');
% Set the title with full LaTeX interpretation
title(t_layout, sprintf('Dual Koopman $y_t^i = C_h(\\kappa_{x_t^i})$ at Time %d ($N_d = %d$)', ...
    target_time, r), 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');

for k = 1:length(ensembles_to_check)
    current_ens = ensembles_to_check(k);
    
    %% =====================================================================
    %% TRUE OPTION 1: THE DUAL/ADJOINT STATE REPRESENTATION
    %% Mathematically evaluates: y_t^i = C_h( f_t^i ) = KMs_r * \kappa_{x_t^i}^\top
    %% =====================================================================
    
    % 1. Calculate the exact empirical snapshot row index for this specific ensemble member
    snapshot_idx = (current_ens - 1) * T_snap + target_time;
    
    % 2. Isolate the exact Dual Koopman State representation from the KEFs matrix
    snapshot_eigenfunctions = KEFs_r(snapshot_idx, :); 
    
    % 3. Project back to physical coordinates via the Output Operator (KMs_r)
    % This computes the RKHS inner products: < h_p, \kappa_{x_{25}^i} >
    spatial_field_vector = real(KMs_r * snapshot_eigenfunctions.'); 
    
    % 4. Reshape back to the true physical [97 x 40] grid geometry
    spatial_field_reshaped = reshape(spatial_field_vector, grid_size);
    
    % --- ACTIVATE NEXT TILE IN LAYOUT ---
    ax = nexttile;
    imagesc(data.y, data.z, spatial_field_reshaped);
    colormap(ax, brighten(redblueTecplot(21), -0.55));
    colorbar;
    axis xy;
    axis equal; % Stop the horizontal/vertical stretching
    xlim([0 2e4]); ylim([0 2e4]);
    
    set(gca, 'FontSize', 11);
    xlabel("y", 'FontSize', 10);
    ylabel("z", 'FontSize', 10);
    title(sprintf('Ensemble #%d', current_ens), 'FontSize', 12, 'FontWeight', 'bold');
    
    %% --- Coordinate Marker Overlay ---
    highlight_row = 25; highlight_col = 26;
    highlight_y = data.y(highlight_col); 
    highlight_z = data.z(highlight_row); 
    hold on;
    plot(highlight_y, highlight_z, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow');
    hold off;
end

set(gcf, 'Renderer', 'painters');



%% Anomonly detetion of ensmebl index 
%% =====================================================================
%% 1. SYSTEM SETUP & ALLOCATION
%% =====================================================================
T_snap = 36;                             % Snapshots per ensemble member
total_snapshots = size(KEFs, 1);         
M = total_snapshots / T_snap;            % Dynamically calculates M=10
grid_points = size(KMs, 1);

% Allocate time-series storage arrays [T_snap x M]
norms_opt1 = zeros(T_snap, M);
norms_opt2 = zeros(T_snap, M);

% Step A: Precompute the full physical true state trajectories
Y_opt1_all = real(KEFs * KMs.'); 

%% =====================================================================
%% 2. CORE MATHEMATICAL CALCULATION (SINGLE PASS LOOP)
%% =====================================================================
for t = 1:T_snap
    % Identify matching frame index 't' across all 10 members
    all_members_at_t = (0:(M-1)) * T_snap + t;
    
    % Compute the Ensemble Mean Fields for step 't'
    mean_phi_t = mean(KEFs(all_members_at_t, :), 1);         % Feature space
    Y_mean_field_t = mean(Y_opt1_all(all_members_at_t, :), 1); % Physical space
    
    for i = 1:M
        row_idx = (i-1)*T_snap + t;
        
        % --- OPTION 1: True Dynamic Deviation ---
        Y_member_t = Y_opt1_all(row_idx, :);
        norms_opt1(t, i) = norm(Y_member_t - Y_mean_field_t);
        
        % --- OPTION 2: Spectral Eigenvalue Propagated Deviation ---
        % Isolate feature space dual amplitude deviation and scale by eigenvalues
        phi_dev_t = KEFs(row_idx, :) - mean_phi_t;
        propagated_amplitude_t = phi_dev_t .* (KEs(:).').^(t-1);
        
        Y_pred_anomaly_t = real(propagated_amplitude_t * KMs.');
        norms_opt2(t, i) = norm(Y_pred_anomaly_t);
    end
end

% Compute total cumulative anomaly tracking scores (Frobenius Equivalent)
deviation_opt1 = sqrt(sum(norms_opt1.^2, 1))'; % Size: M x 1
deviation_opt2 = sqrt(sum(norms_opt2.^2, 1))'; % Size: M x 1

[~, abnormal_idx1] = max(deviation_opt1);
[~, abnormal_idx2] = max(deviation_opt2);


%% =====================================================================
%% VISUALIZATION STYLE A: TIME-EVOLUTION PLUME DIAGRAM (GREY VS RED)
%% =====================================================================
time_steps = 1:T_snap;
figure('Name', 'Temporal Plume Evolution', 'Position', [100, 100, 1100, 420]);

% Left Subplot: Option 1 Plume
subplot(1,2,1); hold on; grid on;
for i = 1:M
    if i ~= abnormal_idx1
        plot(time_steps, norms_opt1(:, i), 'Color', [0.75 0.75 0.75], 'LineWidth', 1.2);
    end
end
plot(time_steps, norms_opt1(:, abnormal_idx1), 'Color', [0.8 0.2 0.2], 'LineWidth', 3);
xlim([1, T_snap]); set(gca, 'FontSize', 11);
xlabel('Time Step (t)'); ylabel('Amplitude Deviation Norm');
title('Opt 1: True Trajectory History');
legend('Normal Cluster', sprintf('Abnormal Member (%d)', abnormal_idx1), 'Location', 'northwest');

% Right Subplot: Option 2 Plume
subplot(1,2,2); hold on; grid on;
for i = 1:M
    if i ~= abnormal_idx2
        plot(time_steps, norms_opt2(:, i), 'Color', [0.75 0.75 0.75], 'LineWidth', 1.2);
    end
end
plot(time_steps, norms_opt2(:, abnormal_idx2), 'Color', [0.9 0.4 0.1], 'LineWidth', 3);
xlim([1, T_snap]); set(gca, 'FontSize', 11);
xlabel('Time Step (t)'); ylabel('Amplitude Deviation Norm');
title('Opt 2: Predicted Spectral History');
legend('Normal Cluster', sprintf('Predicted Outlier (%d)', abnormal_idx2), 'Location', 'northwest');


%% =====================================================================
%% VISUALIZATION STYLE B: CUMULATIVE TRAJECTORY BAR CHARTS
%% =====================================================================
figure('Name', 'Cumulative Integrated Deviation', 'Position', [150, 150, 1000, 420]);

% Left Subplot: Option 1 Bar
subplot(1,2,1);
bar(1:M, deviation_opt1, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none'); hold on;
bar(abnormal_idx1, deviation_opt1(abnormal_idx1), 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none');
grid on; set(gca, 'XTick', 1:M, 'FontSize', 11);
xlabel('Ensemble Index'); ylabel('Cumulative Space-Time Norm');
title('Opt 1: Cumulative True Deviation');

% Right Subplot: Option 2 Bar
subplot(1,2,2);
bar(1:M, deviation_opt2, 'FaceColor', [0.2 0.6 0.5], 'EdgeColor', 'none'); hold on;
bar(abnormal_idx2, deviation_opt2(abnormal_idx2), 'FaceColor', [0.9 0.4 0.1], 'EdgeColor', 'none');
grid on; set(gca, 'XTick', 1:M, 'FontSize', 11);
xlabel('Ensemble Index'); ylabel('Cumulative Space-Time Norm');
title('Opt 2: Cumulative Spectral Deviation');




%% Cehck the low-order model use boxplot--
%% --- Statistical Evaluation of Whole Ensemble Data Matrix ---
% 1. Sort components by dominance
%[~, sort_idx] = sort(abs(KEs), 'descend');
KEs_sorted  = KEs(sort_idx);
KEFs_sorted = KEFs(:, sort_idx);
KMs_sorted  = KMs(:, sort_idx);

r_values = [5, 10, 20, 30, 40, 60, 80, 100, 150]; 
num_r = length(r_values);
total_snapshots = size(Ya_all, 1); % Every row is an ensemble snapshot

% Matrix to store log-residuals for every snapshot across different ranks
log_residual_matrix = zeros(total_snapshots, num_r);

%% 2. One-Shot Matrix Reconstruction for the Whole Ensemble Data
for r_idx = 1:num_r
    r = r_values(r_idx);
    
    % Slice components up to rank r
    KMs_r   = KMs_sorted(:, 1:r);
    KEs_r   = KEs_sorted(1:r);
    KEFs_r  = KEFs_sorted(:, 1:r);
    
    % Reconstruct the WHOLE ensemble data matrix at once
    Ya_pred_r = real(KMs_r * diag(KEs_r) * KEFs_r.');
    
    % 3. Calculate the residual error vector for each individual snapshot row
    for t = 1:total_snapshots
        y_true_snap = Ya_all(t, :);
        y_pred_snap = Ya_pred_r(t, :);
        
        % Normalized error for this specific ensemble row snapshot
        res_norm = norm(y_true_snap - y_pred_snap, 2) / norm(y_true_snap, 2);
        
        % Protect against log10(0) if error becomes exactly zero
        if res_norm < 1e-12, res_norm = 1e-12; end
        
        log_residual_matrix(t, r_idx) = log2(res_norm);
    end
end

%% --- 4. Plotting the Clean Bound Boxplot ---
figure('Position', [100, 100, 720, 440]);

% Generate the boxplot across the whole ensemble data spectrum
boxplot(log_residual_matrix, 'Labels', cellstr(string(r_values)), 'Widths', 0.5);

grid on;
set(gca, 'GridLineStyle', ':', 'LineWidth', 1.2, 'FontSize', 12);
xlabel('Surrogate Model Order (Rank $r$)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Log Normalized Residual $\log_{10}(\|R_{t}^i\|)$', 'Interpreter', 'latex', 'FontSize', 14);
title('Statistical Error Bounds of the Whole Ensemble Dataset $\mathbf{Y}_{\mathrm{all}}$', ...
    'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');

























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




%% Additional tetsing for the KeDMD, K_adj, and K_res 
%%Primal EDMD
% K_edmd = pinv(G1)*A1;
% lam_edmd = eig(K_edmd);
% 
% % Dual operator
% K_adj = A1*pinv(G1);
% lam_dual = eig(K_adj);
% 
% % ResDMD / whitened operator
% lam_res = eig(K_star);
% 
% theta = linspace(0,2*pi,500);
% 
% figure;
% 
% %%-----------------------------
% subplot(1,3,1)
% scatter(real(lam_edmd),imag(lam_edmd),40,'filled');
% hold on
% plot(cos(theta),sin(theta),'k','LineWidth',1.2)
% axis equal
% xlim([-1.2 1.2]); ylim([-1.2 1.2]);
% title('Primal EDMD')
% xlabel('Re(\lambda)')
% ylabel('Im(\lambda)')
% grid on
% 
% %%-----------------------------
% subplot(1,3,2)
% scatter(real(lam_dual),imag(lam_dual),40,'filled');
% hold on
% plot(cos(theta),sin(theta),'k','LineWidth',1.2)
% axis equal
% xlim([-1.2 1.2]); ylim([-1.2 1.2]);
% title('Dual')
% xlabel('Re(\lambda)')
% ylabel('Im(\lambda)')
% grid on
% 
% %%-----------------------------
% subplot(1,3,3)
% scatter(real(lam_res),imag(lam_res),40,'filled');
% hold on
% plot(cos(theta),sin(theta),'k','LineWidth',1.2)
% axis equal
% xlim([-1.2 1.2]); ylim([-1.2 1.2]);
% title('Whitened / ResDMD')
% xlabel('Re(\lambda)')
% ylabel('Im(\lambda)')
% grid on