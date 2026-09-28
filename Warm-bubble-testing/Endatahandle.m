% % Define paths and parameters
baseDir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Dual Koopman-Ensemble-260205\Warm-bubble-testing\EnVorticityData\';
savePathX = fullfile(baseDir, 'EnVfull_Xa.mat');
savePathY = fullfile(baseDir, 'EnVfull_Ya.mat');
%savePath = fullfile(baseDir, 'EnVfull_all260708.mat');

M = 10;          % Number of ensembles
p = 3880;        % Number of spatial points
Tsnap = 51;      % Base snapshot truncation limit
N_sub = Tsnap - 1; % 50 snapshots per ensemble for X and Y

% Preallocate the matrices: 3880 x 500
EnVfull_Xa = zeros(p, M * N_sub);
EnVfull_Ya = zeros(p, M * N_sub);

% Loop through each ensemble file
for i = 0:M-1
    % Construct the file name (Vfull000.mat to Vfull009.mat)
    fileName = sprintf('Vfull%03d.mat', i);
    filePath = fullfile(baseDir, fileName);

    fprintf('Processing: %s\n', fileName);

    % Load the file
    data = load(filePath);

    % Extract the matrix safely
    if isfield(data, 'Vfull')
        current_matrix = data.Vfull;
    else
        varNames = fieldnames(data);
        current_matrix = data.(varNames{1});
    end

    % Extract X (1 to 50) and Y (2 to 51) snapshots for the current ensemble
    X_ensemble = current_matrix(:, 1:N_sub);
    Y_ensemble = current_matrix(:, 2:Tsnap);

    % Calculate column indices for insertion into the master matrices
    startCol = (i * N_sub) + 1;
    endCol   = (i + 1) * N_sub;

    % Insert into master matrices
    EnVfull_Xa(:, startCol:endCol) = X_ensemble;
    EnVfull_Ya(:, startCol:endCol) = Y_ensemble;
end

% Save both matrices separately
save(savePathX, 'EnVfull_Xa', '-v7.3');
save(savePathY, 'EnVfull_Ya', '-v7.3');

fprintf('\nSuccessfully saved:\n1. %s\n2. %s\n', savePathX, savePathY);
%save(savePath, 'EnVfull_all', '-v7.3'); 
fprintf('Successfully saved combined dataset to: %s\n', savePath);


% Ensemable MEan data handle 

% %%Parameters
% p = 3880;
% T_snap = 37;
% M = 10;
% 
% % %1. Reshape the existing 3880 x 500 matrices into 3D arrays
% % Dimensions: [spatial_points, time_snapshots, ensemble_members]
% EnV_Xa_3D = reshape(EnVfull_Xa_36, [p, T_snap, M]);
% EnV_Ya_3D = reshape(EnVfull_Ya_36, [p, T_snap, M]);
% 
% % 2. Calculate the Ensemble Mean across the 3rd dimension (ensemble members)
% % Resulting size will be 3880 x 50
% EnMVfull_Xa_36 = mean(EnV_Xa_3D, 3);
% EnMVfull_Ya_36 = mean(EnV_Ya_3D, 3);
% 
% % 3. Save the resulting mean data
% save(fullfile(baseDir, 'EnMVfull_Xa.mat'), 'EnMVfull_Xa_36', '-v7.3');
% save(fullfile(baseDir, 'EnMVfull_Ya.mat'), 'EnMVfull_Ya_36', '-v7.3');
% 
% fprintf('Ensemble Mean matrices derived successfully.\n');
% fprintf('New dimensions: %d x %d\n', size(EnMVfull_Xa_36, 1), size(EnMVfull_Xa_36, 2));