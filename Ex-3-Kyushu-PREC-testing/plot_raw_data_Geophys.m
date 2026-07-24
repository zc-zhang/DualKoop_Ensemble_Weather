%% plot the raw data_Geophys

%% ============================================================
% LOAD DATA AND INFER DIMENSIONS AUTOMATICALLY
%% ============================================================
%  %data_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Ensemble_MSLP';
% data_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Ensemble_PREC';
% M        = 100;   % number of ensemble members (known from file naming)

%% Load ALL ensemble members
% Xa_ens = cell(M, 1);
% fprintf('Loading %d ensemble members...\n', M);
% for m = 0:M-1
%     %fname        = sprintf('enMSLP%03d.mat', m);
%      fname        = sprintf('enPREC%03d.mat', m);
%     S_load       = load(fullfile(data_dir, fname));
%    % Xraw         = S_load.vectorized_MSLP;        % p × (T+2) unknown size
%     Xraw         = S_load.vectorized_PREC;        % p × (T+2) unknown size
%     Xa_ens{m+1}  = Xraw(:, 2:end);             % p × T: remove boundary cols
% end

%% Load ensemble mean MSLP
% fprintf('Loading mean trajectory (enMSLP100)...\n');
% S_mean    = load(fullfile(data_dir, 'enMSLP100.mat'));
% Xraw_mean = S_mean.vectorized_MSLP;
% Xa_mean   = Xraw_mean(:, 1:end);               % p × T
% 

% %% Load ensemble mean PREC
% fprintf('Loading mean trajectory (enPREC100)...\n');
% S_mean    = load(fullfile(data_dir, 'enPREC100.mat'));
% Xraw_mean = S_mean.vectorized_PREC;
% Xa_mean   = Xraw_mean(:, 2:end);               % p × T
% 
% %% Infer ALL dimensions from loaded data
% p       = size(Xa_ens{1}, 1);    % spatial dimension  = 20860
% T       = size(Xa_ens{1}, 2);    % time snapshots     = 71
% %N_pairs = T - 1;                  % snapshot pairs     = 70
% N_pairs =  T-1;                  % snapshot pairs     = 70
% N_total = M * N_pairs;            % total pairs        = 7100

% fprintf('\n=== Data dimensions ===\n');
% %fprintf('p       = %d  (spatial points = %d x %d grid)\n', p, latN, lonN);
% fprintf('T       = %d  (time snapshots per member)\n',      T);
% fprintf('M       = %d  (ensemble members)\n',               M);
% fprintf('N_pairs = %d  (snapshot pairs per member)\n',      N_pairs);
% fprintf('N_total = %d  (total snapshot pairs, M*N_pairs)\n',N_total);
% fprintf('Xa_mean size: %d x %d\n', size(Xa_mean,1), size(Xa_mean,2));
% 
% %% Verify consistency across members
% % Check all members have same dimensions as member 1
% for m = 1:M
%     assert(size(Xa_ens{m},1)==p, 'Member %d: wrong spatial dim',m);
%     assert(size(Xa_ens{m},2)==T, 'Member %d: wrong time dim',m);
% end
% assert(size(Xa_mean,1)==p, 'Mean: wrong spatial dim');
% assert(size(Xa_mean,2)==T, 'Mean: wrong time dim');
% fprintf('All %d members verified: consistent dimensions.\n', M);
% 
% %% ============================================================
% % STEP 1: STACK ALL ENSEMBLES (Scenario A, fat-wide)
% % Build Xa_all, Ya_all from DEVIATIONS from mean
% % Xa_all: p × N_total = 20860 × 7000
% % Ya_all: p × N_total = 20860 × 7000
% %% ============================================================
% 
% % Pre-allocate using inferred dimensions (faster than growing array)
% Xa_all = zeros(p, N_total);
% Ya_all = zeros(p, N_total);
% 
% % In your loading loop:
% col = 1;
% for m = 1:M
%     % Use the raw snapshots directly
%     Xa_all(:, col:col+N_pairs-1) = Xa_ens{m}(:, 2:end-1);  % Raw p × N_pairs
%     Ya_all(:, col:col+N_pairs-1) = Xa_ens{m}(:, 3:end);    % Raw p × N_pairs
%     col = col + N_pairs;
% end
% % 
% % % % Similarly for initial conditions
% X0_all = zeros(p, M);   % pre-allocate with inferred p, M
% for m = 1:M
%     X0_all(:,m) = Xa_ens{m}(:,1); % Raw initial state
% end


%=================Create New Ensemble Data set==================
% % OPTION C: Deviation data (X^i -\barX)================
% col = 1;
% for m = 1:M
%     % Deviation from ensemble mean
%     dev_m = Xa_ens{m} - Xa_mean;               % p × T
% 
%     % Snapshot pairs from deviation trajectory
%     Xa_all(:, col:col+N_pairs-1) = dev_m(:, 1:end-1);  % p × N_pairs
%     Ya_all(:, col:col+N_pairs-1) = dev_m(:, 2:end);    % p × N_pairs
% 
%     col = col + N_pairs;
% end

% fprintf('\nXa_all: %d × %d (p × N_total)\n', size(Xa_all,1), size(Xa_all,2));
% fprintf('Ya_all: %d × %d (p × N_total)\n', size(Ya_all,1), size(Ya_all,2));
% 
% ===========Initial Deviation=========
% %========
% X0_all = zeros(p, M);   % pre-allocate with inferred p, M
% for m = 1:M
%     dev_m       = Xa_ens{m} - Xa_mean;
%     X0_all(:,m) = dev_m(:,1);   % deviation at first time snapshot
% end
% 
% fprintf('X0_all: %d × %d (p × M, initial conditions)\n', ...
%     size(X0_all,1), size(X0_all,2));
%============================================================================

% ============================================================
%STEP 3: SAVE (using inferred variable names with dimensions)
% ============================================================
% output_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Dual Koopman-Ensemble-260205\data';
% if ~exist(output_dir,'dir'); mkdir(output_dir); end
% 
% output_file = fullfile(output_dir, 'EnPREC_all_col_raw.mat');
% fprintf('\nSaving to: %s\n', output_file);
% save(output_file, ...
%     'Xa_all', 'Ya_all', ...       % fat-wide Scenario A matrices
%     'X0_all', ...                  % initial conditions
%     'Xa_mean', ...                 % ensemble mean trajectory
%     'Xa_ens', ...                  % all individual members (cell)
%     'p', 'T', 'M', 'N_pairs', 'N_total', ...  % dimensions
%     '-v7.3');
% fprintf('Saved successfully.\n');
% 
% % Print final summary
% fprintf('\n=== Final data summary ===\n');
% fprintf('Xa_all:  %d × %d  (deviation snapshot pairs, Scenario A)\n',...
%     size(Xa_all,1), size(Xa_all,2));
% fprintf('Ya_all:  %d × %d  (one-step-ahead deviations)\n',...
%     size(Ya_all,1), size(Ya_all,2));
% fprintf('X0_all:  %d × %d  (initial condition deviations)\n',...
%     size(X0_all,1), size(X0_all,2));
% fprintf('Xa_mean: %d × %d  (ensemble mean trajectory)\n',...
%     size(Xa_mean,1), size(Xa_mean,2));
% fprintf('G1 will be: %d × %d\n', N_total, N_total);
% fprintf('Memory of G1: %.1f MB\n', N_total^2*8/1e6);
% 




%%%%% Diretcly lodad data (datahandle)
%load('weatherDat2021AUG_ensemble0.mat');
%load('EnPREC_all_col_raw.mat');
%load('EnMSLP_all_col_raw.mat');



%% collect 48 hrs testing data
%% 1. Parameters (Updated)
% %data_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Ensemble_PREC';
% data_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Ensemble_MSLP';
% M        = 100;   
% N_pairs  = 50;    
% p        = 20860; 
% 
% %% 2. Loading and Stacking
% Xa_all = zeros(p, M * N_pairs);
% Ya_all = zeros(p, M * N_pairs);
% X0_all = zeros(p, M);
% 
% % Pre-allocate the Mean
% Xa_mean_processed = zeros(p, N_pairs);
% 
% % Process Ensemble Members
% for m = 1:M
%     %fname = sprintf('enPREC%03d.mat', m-1); 
%     fname = sprintf('enMSLP%03d.mat', m-1); 
%     S_load = load(fullfile(data_dir, fname));
%     %Xraw = S_load.vectorized_PREC(:, 2:end); 
%     Xraw = S_load.vectorized_MSLP(:, 2:end); 
% 
%     start_idx = (m-1) * N_pairs + 1;
%     end_idx   = m * N_pairs;
% 
%     Xa_all(:, start_idx:end_idx) = Xraw(:, 1:N_pairs);
%     Ya_all(:, start_idx:end_idx) = Xraw(:, 2:N_pairs+1);
%     X0_all(:, m) = Xraw(:, 1);
% end
% 
% %% 3. Process the Ensemble Mean (enPREC100)
% S_mean = load(fullfile(data_dir, 'enMSLP100.mat'));
% %S_mean = load(fullfile(data_dir, 'enPREC100.mat'));
% % Assuming the mean file has the same structure as the others
% %Xraw_mean = S_mean.vectorized_PREC(:, 2:end); 
% Xraw_mean = S_mean.vectorized_MSLP(:, 2:end); 
% Xa_mean = Xraw_mean(:, 1:N_pairs);
% 
% %% 4. Saving
% % Save the mean alongside your other matrices
% %save('EnPREC_all_col_raw_48hrs.mat', 'Xa_all', 'Ya_all', 'X0_all', 'Xa_mean', '-v7.3');
% save('EnMSLP_all_col_raw_48hrs.mat', 'Xa_all', 'Ya_all', 'X0_all', 'Xa_mean', '-v7.3');
% fprintf('Successfully processed all ensembles and the mean. Data saved.\n');



%%   =======72 hors
% %% 1. Parameters
% data_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Ensemble_PREC';
% M        = 100;   
% N_pairs  = 70;    
% p        = 20860; 
% 
% %% 2. Pre-allocation
% Xa_all = zeros(p, M * N_pairs);
% Ya_all = zeros(p, M * N_pairs);
% X0_all = zeros(p, M);
% Xa_ens = cell(M, 1);
% 
% fprintf('Processing PREC ensembles...\n');
% 
% %% 3. Main Loading Loop
% for m = 1:M
%     fname = sprintf('enPREC%03d.mat', m-1); 
%     S_load = load(fullfile(data_dir, fname));
%     Xraw = S_load.vectorized_PREC; 
% 
%     % The "Standard" indexing:
%     % Xa uses columns 2 through 50 (Total 49)
%     % Ya uses columns 3 through 51 (Total 49)
%     Xa_ens{m} = Xraw(:, 3 : 3 + N_pairs - 1);
% 
%     start_idx = (m-1) * N_pairs + 1;
%     end_idx   = m * N_pairs;
% 
%     Xa_all(:, start_idx:end_idx) = Xraw(:, 3 : 3 + N_pairs - 1);
%     Ya_all(:, start_idx:end_idx) = Xraw(:, 4 : 4 + N_pairs - 1);
% 
%     % Set X0 to the first column of the data we are actually using (Column 2)
%     X0_all(:, m) = Xraw(:, 3); 
% end
% 
% %% 4. Process the Ensemble Mean
% Xa_mean = zeros(p, N_pairs);
% for m = 1:M
%     Xa_mean = Xa_mean + Xa_ens{m};
% end
% Xa_mean = Xa_mean / M;
% 
% % 5. Saving
% % save('EnPREC_all_col_raw_72hrs.mat', 'Xa_all', 'Ya_all', 'X0_all', 'Xa_mean', 'Xa_ens', '-v7.3');
% % fprintf('Successfully saved PREC data to EnPREC_all_col_raw_72hrs.mat\n');
% 
% 
% output_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Dual Koopman-Ensemble-260205\data';
% if ~exist(output_dir,'dir'); mkdir(output_dir); end
% 
% output_file = fullfile(output_dir, 'EnPREC_all_col_raw_36hrs.mat');
% fprintf('\nSaving to: %s\n', output_file);
% save(output_file, ...
%     'Xa_all', 'Ya_all', ...       % fat-wide Scenario A matrices
%     'X0_all', ...                  % initial conditions
%     'Xa_mean', ...                 % ensemble mean trajectory
%     'Xa_ens', ...                  % all individual members (cell)
%     'p', 'T', 'M', 'N_pairs', 'N_total', ...  % dimensions
%     '-v7.3');
% fprintf('Saved successfully.\n');



%%  36 hrs 
% %% 1. Parameters
% data_dir   = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Ensemble_PREC';
% output_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Dual Koopman-Ensemble-260205\data';
% M          = 100;   
% N_pairs    = 36;    
% p          = 20860; 
% T          = 36;              % Total duration in hours
% N_total    = M * N_pairs;     % Total number of snapshots
% 
% %% 2. Pre-allocation
% Xa_all  = zeros(p, N_total);
% Ya_all  = zeros(p, N_total);
% X0_all  = zeros(p, M);
% Xa_ens  = cell(M, 1);
% 
% fprintf('Processing PREC ensembles for 36-hour trajectory...\n');
% 
% %% 3. Main Loading Loop
% for m = 1:M
%     fname = sprintf('enPREC%03d.mat', m-1); 
%     S_load = load(fullfile(data_dir, fname));
%     Xraw = S_load.vectorized_PREC; 
% 
%     % Define the sliding window columns
%     % Xa: Column 3 to 3+34-1 = 36
%     % Ya: Column 4 to 4+34-1 = 37
%     idx_Xa = 3 : (3 + N_pairs - 1);
%     idx_Ya = 4 : (4 + N_pairs - 1);
% 
%     % Store members in cell array
%     Xa_ens{m} = Xraw(:, idx_Xa);
% 
%     % Map to global fat-wide matrices
%     start_idx = (m-1) * N_pairs + 1;
%     end_idx   = m * N_pairs;
% 
%     Xa_all(:, start_idx:end_idx) = Xraw(:, idx_Xa);
%     Ya_all(:, start_idx:end_idx) = Xraw(:, idx_Ya);
% 
%     % Initial condition (start of each ensemble member)
%     X0_all(:, m) = Xraw(:, 3); 
% end
% 
% %% 4. Process the Ensemble Mean
% % Calculate mean along the 3rd dimension of the cell contents
% % Xa_mean = mean(cat(3, Xa_ens{:}), 3);
% % %% 4. Process the Ensemble Mean
% Xa_mean = zeros(p, N_pairs);
% for m = 1:M
%     Xa_mean = Xa_mean + Xa_ens{m};
% end
% Xa_mean = Xa_mean / M;
% 
% %% 5. Saving the Data
% if ~exist(output_dir,'dir'); mkdir(output_dir); end
% output_file = fullfile(output_dir, 'EnPREC_all_col_raw_36hrs.mat');
% 
% fprintf('\nSaving to: %s\n', output_file);
% save(output_file, ...
%     'Xa_all', 'Ya_all', ...       % Koopman Data Matrices
%     'X0_all', ...                  % Initial Conditions
%     'Xa_mean', ...                 % Ensemble Mean
%     'Xa_ens', ...                  % Ensemble Members
%     'p', 'T', 'M', 'N_pairs', 'N_total', ...  % Metadata
%     '-v7.3');
% 
% fprintf('Data generation and saving complete.\n');
%% Geophysical grid information
lat_sub = weatherDat2021AUG_ensemble0.dat_lat;  % 140×149
lon_sub = weatherDat2021AUG_ensemble0.dat_lon;  % 140×149
latN = 140;
lonN = 149;
shapefile_path = 'D:\Susuki Lab\Testing_Code\data-weather\Data_250401\ne_10m_coastline\ne_10m_coastline.shp';
S = shaperead(shapefile_path);


%% 
% 2. Calculate index for Kumamoto (32.803N, 130.708E)
targetLat = 32.803; targetLon = 130.708;
dist = sqrt((lat_sub - targetLat).^2 + (lon_sub - targetLon).^2);
[minDist, Locatidx] = min(dist(:));
[r, c] = ind2sub(size(lat_sub), Locatidx);
fprintf('Kumamoto is at Index: %d (Row: %d, Col: %d)\n', Locatidx, r, c);


%===============
% 3. Calculate mean
field_mean_map = reshape(sum(Xa_mean, 2), [latN, lonN]);

% 4. Plotting
figure('Color','w','Position',[100 100 800 600]);
pcolor(lon_sub, lat_sub, field_mean_map);
shading interp; 
axis equal tight;
set(gca,'YDir','normal');

% Colormap
colormap(brighten(redblueTecplot(21), -0.55));
colorbar;
hold on;

% Add Coastlines
for kk = 1:length(S)
    plot(S(kk).X, S(kk).Y, 'k-', 'LineWidth', 1.5);
end

%% Add Kumamoto Star Marker
plot(lon_sub(r, c), lat_sub(r, c), 'p', ...
    'MarkerSize', 14, ...
    'MarkerFaceColor', 'g', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1);

hold off;

% 5. Axis Formatting
set(gca, 'TickLabelInterpreter', 'tex', 'YTickLabelRotation', 0);
xtickformat('%g^{\\circ}E');
ytickformat('%g^{\\circ}N');

title('Ensemble mean (all members)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('Longitude', 'FontSize', 11);
ylabel('Latitude', 'FontSize', 11);

xlim([min(lon_sub(:)), max(lon_sub(:))]);
ylim([min(lat_sub(:)), max(lat_sub(:))]);



%% ============raw data plot (timne series）

% 1. Create the time vector (72 hours starting Aug 10, 12:00, 2021)
startTime = datetime(2021, 8, 10, 12, 0, 0);
timeVec = startTime + hours(0:35);   % PREC 71  %MSLP 72

% 2. Calculate the mean over the ensemble
dataToPlot = mean(Xa_mean, 1);
dataToPlot_Locatindx = Xa_mean(Locatidx,:);

% 3. Plotting
figure('Color', 'w', 'Position', [100 100 900 450]);
bar(timeVec, dataToPlot, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
hold on;
bar(timeVec, dataToPlot_Locatindx, 'FaceColor', [0.1 0.5 0.4], 'EdgeColor', 'none');


% Added Legend
legend({'Ensemble Mean of Whole Kyushu Field', 'Ensemble Mean of Kumamoto'}, ...
    'Location', 'best', 'FontSize', 14);

% 4. Customizing Axis Ticks
ax = gca;
tickTimes = startTime:hours(6):(startTime + hours(36));
ax.XTick = tickTimes;

% 5. Build two-line tick labels for R2020b
tickLabels = cell(numel(tickTimes), 1);
for i = 1:numel(tickTimes)
    t = tickTimes(i);
    timeStr = datestr(t, 'HH:MM');
    dateStr = datestr(t, 'mm/dd');
    tickLabels{i} = [timeStr, ' ', dateStr];
end

% Set as categorical to force string labels on datetime axis
ax.XAxis.TickValues = tickTimes;
set(ax, 'XTickLabel', tickLabels);
set(ax, 'XTickLabelRotation', 20);

% 6. Aesthetics
grid on;
xlabel('Date/Time (UTC)', 'FontSize', 20);
ylabel('Ensemble Mean', 'FontSize', 20);
title('Ensemble Mean PREC Time Series (Aug 10 - Aug 13, 2021)', 'FontSize', 20);
set(gca, 'FontSize', 20);





%%==================
% PLOT : Ensemble mean MSLP
% Xa_mean is already loaded: p × T = 20860 × 71
% Just pick which time step to show, or show time-averaged mean
%% ============================================================
figure('Color','w','Position',[100 100 800 600]);

%% Option A: mean over all time steps (one map)
%field_mean_all = reshape(mean(Xa_mean,2), [latN, lonN]);   % p×1 -> latN×lonN

field_mean_map = reshape(sum(Xa_mean, 2), [latN, lonN]);   % sum over 72 steps
%% Option B: show at a specific time step k
 %field_mean_all = reshape(Xa_mean(:, k_show+1), [latN, lonN]);

%pcolor(lon_sub, lat_sub, field_mean_all);
pcolor(lon_sub, lat_sub, field_mean_map);
shading interp; axis equal tight;

set(gca,'YDir','normal');

% % Axis Handling
% set(gca,'YDir','normal');
% set(gca, 'TickLabelInterpreter', 'tex', ...
%          'YTickLabelRotation', 90); % Ensures Y-axis labels are horizontal
%colormap(jet); 
colormap(brighten(redblueTecplot(21),-0.55));
colorbar;

hold on;
for kk=1:length(S)
    plot(S(kk).X, S(kk).Y, 'k-', 'LineWidth', 1.5);
end


%% Add Kumamoto Star Marker
plot(lon_sub(r, c), lat_sub(r, c), 'p', ...
    'MarkerSize', 14, ...
    'MarkerFaceColor', 'g', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1);

hold off;

title('Ensemble mean (all members)',...
    'FontSize',13,'FontWeight','bold');
 xlabel('Longitude (^\circE)');
 ylabel('Latitude (^\circN)');
xlim([min(lon_sub(:)), max(lon_sub(:))]);
ylim([min(lat_sub(:)), max(lat_sub(:))]);



%%      =================Alternative P>LOT
figure('Color','w','Position',[100 100 800 600]);

% Data Preparation (As per your provided logic)
field_mean_map = reshape(sum(Xa_mean, 2), [latN, lonN]);

% Plotting
pcolor(lon_sub, lat_sub, field_mean_map);
shading interp; 
axis equal tight;

% Axis Handling
set(gca,'YDir','normal');
set(gca, 'TickLabelInterpreter', 'tex', ...
         'YTickLabelRotation', 90); % Ensures Y-axis labels are horizontal

% Customizing Tick Labels to include degree and cardinal directions
xtickformat('%g^{\\circ}E');
ytickformat('%g^{\\circ}N');

% Colormap and aesthetics
colormap(brighten(redblueTecplot(21), -0.55));
colorbar;
hold on;

% Plotting the lines (S)
for kk = 1:length(S)
    plot(S(kk).X, S(kk).Y, 'k-', 'LineWidth', 1.5);
end


%% Add Kumamoto Star Marker
plot(lon_sub(r, c), lat_sub(r, c), 'p', ...
    'MarkerSize', 14, ...
    'MarkerFaceColor', 'g', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1);

hold off;

%hold off;

% Titles and Axis Labels
title('Ensemble mean (all members)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('Longitude', 'FontSize', 11);
ylabel('Latitude', 'FontSize', 11);

% Final bounds
xlim([min(lon_sub(:)), max(lon_sub(:))]);
ylim([min(lat_sub(:)), max(lat_sub(:))]);


%% %=====================================================================
%% Chose some one ensmeble to see the rawaensmebel (not ensmeble mean data )

figure('Color','w','Position',[100 100 800 600]);

% Data Preparation (As per your provided logic)
% load the one of PREC/MSLP ensemble data 
Xa_raw_ens1=vectorized_PREC(:,3:39-1);  % set the time
Xa_dev =Xa_raw_ens1 -Xa_mean;
%field_raw_map = reshape(sum(Xa_raw_ens1, 2), [latN, lonN]);
field_dev_map = reshape(sum(Xa_dev, 2), [latN, lonN]);
% Plotting
%pcolor(lon_sub, lat_sub, field_raw_map);
pcolor(lon_sub, lat_sub, field_dev_map);
shading interp; 
axis equal tight;

% Axis Handling
set(gca,'YDir','normal');
set(gca, 'TickLabelInterpreter', 'tex', ...
         'YTickLabelRotation', 90); % Ensures Y-axis labels are horizontal

% Customizing Tick Labels to include degree and cardinal directions
xtickformat('%g^{\\circ}E');
ytickformat('%g^{\\circ}N');

% Colormap and aesthetics
colormap(brighten(redblueTecplot(21), -0.55));
colorbar;
hold on;

% Plotting the lines (S)
for kk = 1:length(S)
    plot(S(kk).X, S(kk).Y, 'k-', 'LineWidth', 1.5);
end


%% Add Kumamoto Star Marker
plot(lon_sub(r, c), lat_sub(r, c), 'p', ...
    'MarkerSize', 14, ...
    'MarkerFaceColor', 'g', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1);

hold off;

%hold off;

% Titles and Axis Labels
%title('Total Raw PREC (Ensemble 1, 36 hrs/mm)', 'FontSize', 13, 'FontWeight', 'bold');
title('Total Deviation PREC (Ensemble 1, 36 hrs/mm)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('Longitude', 'FontSize', 11);
ylabel('Latitude', 'FontSize', 11);

% Final bounds
xlim([min(lon_sub(:)), max(lon_sub(:))]);
ylim([min(lat_sub(:)), max(lat_sub(:))]);





%% Bar plot of raw 
% 1. Create the time vector (72 hours starting Aug 10, 12:00, 2021)
startTime = datetime(2021, 8, 10, 12, 0, 0);
timeVec = startTime + hours(0:35);   % PREC 71  %MSLP 72

% 2. Calculate the mean over the ensemble
dataToPlot_raw = mean(Xa_raw_ens1, 1);
dataToPlot_Locatindx_raw = Xa_raw_ens1(Locatidx,:);

% 3. Plotting
figure('Color', 'w', 'Position', [100 100 900 450]);
bar(timeVec, dataToPlot_raw, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
hold on;
bar(timeVec, dataToPlot_Locatindx_raw, 'FaceColor', [0.1 0.5 0.4], 'EdgeColor', 'none');


% Added Legend
legend({'Raw Ensemble one of Whole Kyushu Field', 'Raw PREC of Kumamoto'}, ...
    'Location', 'best', 'FontSize', 14);

% 4. Customizing Axis Ticks
ax = gca;
tickTimes = startTime:hours(6):(startTime + hours(36));
ax.XTick = tickTimes;

% 5. Build two-line tick labels for R2020b
tickLabels = cell(numel(tickTimes), 1);
for i = 1:numel(tickTimes)
    t = tickTimes(i);
    timeStr = datestr(t, 'HH:MM');
    dateStr = datestr(t, 'mm/dd');
    tickLabels{i} = [timeStr, ' ', dateStr];
end

% Set as categorical to force string labels on datetime axis
ax.XAxis.TickValues = tickTimes;
set(ax, 'XTickLabel', tickLabels);
set(ax, 'XTickLabelRotation', 20);

% 6. Aesthetics
grid on;
xlabel('Date/Time (UTC)', 'FontSize', 20);
ylabel('Mean', 'FontSize', 20);
title('Ensemble Index 1 of Raw PREC Time Series (Aug. 10 - Aug. 13, 2021)', 'FontSize', 20);
set(gca, 'FontSize', 20);


%% Deviation of esemble mean 
% Chose some one ensmeble to see the rawaensmebel (not ensmeble mean data )

figure('Color','w','Position',[100 100 800 600]);

% Data Preparation (As per your provided logic)
% load the one of PREC/MSLP ensemble data 
%Xa_raw_ens1=vectorized_PREC(:,3:39-1);  % set the time
Xa_dev =Xa_raw_ens1 -Xa_mean;
%field_raw_map = reshape(sum(Xa_raw_ens1, 2), [latN, lonN]);
field_dev_map = reshape(sum(Xa_dev, 2), [latN, lonN]);
%%Plotting
%%pcolor(lon_sub, lat_sub, field_raw_map);
pcolor(lon_sub, lat_sub, field_dev_map);
shading interp; 
axis equal tight;

% Axis Handling
set(gca,'YDir','normal');
set(gca, 'TickLabelInterpreter', 'tex', ...
         'YTickLabelRotation', 90); % Ensures Y-axis labels are horizontal

% Customizing Tick Labels to include degree and cardinal directions
xtickformat('%g^{\\circ}E');
ytickformat('%g^{\\circ}N');

% Colormap and aesthetics
%colormap(brighten(redblueTecplot(21), -0.55));
colormap(jet);
colorbar;
hold on;

% Plotting the lines (S)
for kk = 1:length(S)
    plot(S(kk).X, S(kk).Y, 'k-', 'LineWidth', 1.5);
end


%% Add Kumamoto Star Marker
plot(lon_sub(r, c), lat_sub(r, c), 'p', ...
    'MarkerSize', 14, ...
    'MarkerFaceColor', 'g', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1);

hold off;

%hold off;

% Titles and Axis Labels
%title('Total Raw PREC (Ensemble 1, 36 hrs/mm)', 'FontSize', 13, 'FontWeight', 'bold');
title('Total Deviation PREC (Ensemble 1, 36 hrs/mm)', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('Longitude', 'FontSize', 11);
ylabel('Latitude', 'FontSize', 11);

% Final bounds
xlim([min(lon_sub(:)), max(lon_sub(:))]);
ylim([min(lat_sub(:)), max(lat_sub(:))]);


