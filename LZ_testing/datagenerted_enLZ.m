% LZ mode data geenrating ==

%% -------------------- STEP 1: Lorenz Ensemble Generation --------------------
clear; clc;

% SWITCH CONFIGURATION HERE:
% 'A' = Climate / Data Assimilation (Tight Local Cluster, tracks Ensemble Mean)
% 'B' = Pure Koopman Spectral (Global Distributed Spread across both wings)

% trajectories are generated every time you run the script.
rng(42);

RUN_MODE = 'B'; 

SIGMA = 10; RHO = 28; BETA = 8/3;   % Classic chaotic Lorenz parameters
ODEFUN = @(t,y) [SIGMA*(y(2)-y(1));
                 y(1).*(RHO-y(3))-y(2);
                 y(1).*y(2)-BETA*y(3)];
             
options = odeset('RelTol',1e-10,'AbsTol',1e-11);
delta_t = 0.05;   

%% ---- SET CONDITIONAL PARAMETERS BASED ON MODE ----
if strcmp(RUN_MODE, 'A')
    fprintf('=== RUNNING TEST A: Climate / Data Assimilation (Local Tracking) ===\n');
    M       = 50;     % Preferred ensemble size for localized cloud tracking
    T_snap  = 65;    % Short snapshots (ideal for tracking local evolution)
    
    % Initialize all members in a tight ball on ONE WING center
    base_IC = [-10.0; -15.0; 25.0]; 
    X0_raw  = 0.015 * (rand(3, M) - 0.5) + base_IC;
    
else
    fprintf('=== RUNNING TEST B: Pure Koopman Spectral (Global Diversity) ===\n');
    M       = 20;     % Fewer members, but longer tracks
    T_snap  = 1000;   % Much longer tracks to evaluate global operator structure
    
    % Spread ICs across a MUCH larger region anchoring both wings
    C_plus  = [ sqrt(BETA*(RHO-1));  sqrt(BETA*(RHO-1)); RHO-1];
    C_minus = [-sqrt(BETA*(RHO-1)); -sqrt(BETA*(RHO-1)); RHO-1];
    half    = floor(M/2);
    X0_raw  = [ 3*(rand(3,half)  -0.5) + C_plus, ...
                3*(rand(3,M-half)-0.5) + C_minus ];
end

%% ---- SIMULATION LOOP ----
Xa_ens = cell(M,1); Ya_ens = cell(M,1); X0_all = zeros(3,M);

for m = 1:M
    if strcmp(RUN_MODE, 'A')
        % Fixed burn-in keeps the tight localized cloud intact at evaluation start
        T_burn_eff = 5; 
    else
        % Randomized burn-in decorrelates phases globally in time and space
        T_burn_eff = 5 + 15*rand; 
    end
    
    [~,Y_burn] = ode45(ODEFUN, [0 T_burn_eff], X0_raw(:,m), options);
    X0 = Y_burn(end,:).';
    X0_all(:,m) = X0;
    
    t_span = 0:delta_t:(T_snap*delta_t);
    [~,Y_orbit] = ode45(ODEFUN, t_span, X0, options);
    
    Xa_ens{m} = Y_orbit(1:end-1,:).'; 
    Ya_ens{m} = Y_orbit(2:end,:).';   
end

% Pool ALL members for mapping
Xa_all = cell2mat(Xa_ens');
Ya_all = cell2mat(Ya_ens');

%% ---- OPTIONAL FEATURE THINNING (Only for Mode B if Gram matrix is too large) ----
if strcmp(RUN_MODE, 'B')
    % k = 4;
    % Xa_all = Xa_all(:,1:k:end);
    % Ya_all = Ya_all(:,1:k:end);
end

[dim, Total_Snapshots] = size(Xa_all);
fprintf('Done. Total snapshots pooled: %d\n\n', Total_Snapshots);


%% ---- STEP 2: Save Data for Koopman Pipeline ----
% Define a dynamic filename based on the selected mode
filename = sprintf('LZ_ensemble_data_Test_%s.mat', RUN_MODE);

% Save the vital matrices, initial conditions, and parameters
save(filename, 'Xa_all', 'Ya_all', 'X0_all', 'RUN_MODE', 'delta_t', 'M', 'T_snap');

fprintf('Successfully saved %d snapshots to file: %s\n', Total_Snapshots, filename);