% ensembledata_reshape
% EnAug21DataReshape.m

%% Parameters
%data_dir = 'D:\...\Ensemble_PREC';
%  %data_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Ensemble_MSLP';
data_dir = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Ensemble_PREC';
M        = 100;       % ensemble members
N_pairs  = 70;        % snapshot pairs per member
p        = 20860;     % state dimension
lag      = 1;         % time-step lag between Xa and Ya

%% Pre-allocate
Xa_all = zeros(p, M * N_pairs);
Ya_all = zeros(p, M * N_pairs);
X0_all = zeros(p, M);
Xa_ens = cell(M, 1);

%% Load and stack
for m = 1:M
    fname = sprintf('enPREC%03d.mat', m-1);   % 000..099
    S     = load(fullfile(data_dir, fname));
    Xraw  = S.vectorized_PREC;                % p × N_total

    col_start = 3;                            % first column to use
    xa_cols = col_start : col_start + N_pairs - 1;        % 3..72
    ya_cols = col_start + lag : col_start + lag + N_pairs - 1; % 4..73

    Xa_block = Xraw(:, xa_cols);   % p × N_pairs
    Ya_block = Xraw(:, ya_cols);   % p × N_pairs

    idx = (m-1)*N_pairs + 1 : m*N_pairs;
    Xa_all(:, idx) = Xa_block;
    Ya_all(:, idx) = Ya_block;

    X0_all(:, m)   = Xraw(:, col_start);   % initial condition of this member
    Xa_ens{m}      = Xa_block;
end

%% Ensemble mean trajectory  (p × N_pairs)
Xa_mean = mean(cat(3, Xa_ens{:}), 3);

%% Save
output_dir  = 'D:\Susuki Lab\Testing_Code\data-weather\Data250525\Dual Koopman-Ensemble-260205\data';
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

N_total     = N_pairs + lag;          % total time steps used per member

output_file = fullfile(output_dir, 'EnPREC_all_col_raw_36hrs.mat');
fprintf('\nSaving to: %s\n', output_file);

save(output_file, ...
    'Xa_all',   ...    % p × (M*N_pairs)  — all Xa snapshot pairs
    'Ya_all',   ...    % p × (M*N_pairs)  — all Ya snapshot pairs
    'X0_all',   ...    % p × M            — initial conditions per member
    'Xa_mean',  ...    % p × N_pairs      — ensemble mean trajectory
    'Xa_ens',   ...    % M×1 cell, each p × N_pairs
    'p',        ...    % state dimension = 20860
    'M',        ...    % ensemble members = 100
    'N_pairs',  ...    % snapshot pairs per member = 70
    'N_total',  ...    % time steps per member used = 71
    'lag',      ...    % time lag between Xa and Ya = 1
    '-v7.3');

fprintf('Saved successfully.\n');