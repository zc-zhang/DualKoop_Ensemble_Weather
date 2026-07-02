%% -------------------- BOXPLOT: Reconstruction Accuracy vs N_dict --------------------
% Assumes Xa_all, Ya_all already exist in the workspace (pooled Lorenz snapshots),
% and that kernel_ResDMD is on the path (your local version returning
% [G, K, L, PX, PY, PSI_x, PSI_y, PSI_y2, G1, A1, kernel_f]).

%clc; %clear all;
N_dict_list = [10 20 50 100 150 200 300 500];
R_reps      = 15;              % repetitions per N_dict value
frac_sub    = 0.85;            % fraction of pooled snapshots used per repetition

METHOD = 'bootstrap';          % 'bootstrap' (fast) or 'reseed' (regenerate Lorenz ensemble)

% ---- Generalized Lorenz model config (only used when METHOD = 'reseed') ----
CASE.dim     = 6;              % 3, 5, 6, 8, 9, or 11
CASE.SCALE   = 1/10;           % coupling scale (classic 3D case often uses SCALE=1)
CASE.SIGMA   = 10;
CASE.BETA    = 8/3;
CASE.RHO     = 40;             % classic 3D chaos uses RHO~28; higher-dim cases often use RHO~40+
CASE.M       = 50;             % ensemble members per repetition
CASE.delta_t = 0.05;
CASE.T_snap  = 40;
CASE.T_burn  = 5.0;

Total_Snapshots = size(Xa_all, 2);
errs = nan(R_reps, numel(N_dict_list));   % rows = repetitions, cols = N_dict values

fprintf('Running %d repetitions x %d N_dict values (method = %s)...\n', ...
        R_reps, numel(N_dict_list), METHOD);

for k = 1:numel(N_dict_list)
    N_dict = N_dict_list(k);
    for r = 1:R_reps

        switch METHOD
            case 'bootstrap'
                % Resample a random subset of the pooled snapshot pairs to
                % build the dictionary from; test/reconstruct on that same
                % subset for a like-for-like error each repetition.
                rng(1000*k + r);   % reproducible but distinct per (N_dict, rep)
                n_sub = round(frac_sub * Total_Snapshots);
                idx = randperm(Total_Snapshots, n_sub);
                Xa_use = Xa_all(:, idx);
                Ya_use = Ya_all(:, idx);

            case 'reseed'
                % Regenerate a fresh Lorenz ensemble with a new seed, then
                % re-pool. Requires the ensemble-generation code from the
                % top of the original script wrapped as a local function
                % (see generate_lorenz_pool below) -- slower but tests true
                % statistical variability across independent trajectories.
                rng(1000*k + r);
                [Xa_use, Ya_use] = generate_lorenz_pool(CASE);
        end

        try
            [G, K_star, L, PX, PY, PSI_x, PSI_y, PSI_y2, G1, A1, kernel_f] = ...
                kernel_ResDMD(Xa_use, Ya_use, 'type', 'Gaussian', 'N', N_dict);

            % No slicing needed -- G, K_star, PX are already N_dict-sized
            [V_coeff, Lambda_mat] = eig(K_star, G);   % G ~= eye(N_dict)
            KEFs = PX * V_coeff;

            KMs     = Xa_use * pinv(KEFs).';
            Ya_pred = KMs * Lambda_mat * KEFs.';

            errs(r, k) = norm(Ya_use - Ya_pred, 'fro') / norm(Ya_use, 'fro') * 100;
        catch ME
            fprintf('  [N_dict=%d, rep=%d] failed: %s\n', N_dict, r, ME.message);
            errs(r, k) = NaN;
        end
    end
    fprintf('N_dict = %4d  ->  median error = %.3f%%  (range %.3f - %.3f%%)\n', ...
        N_dict, median(errs(:,k), 'omitnan'), min(errs(:,k)), max(errs(:,k)));
end

%% -------------------- PLOT --------------------
figure('Position', [100, 100, 850, 480]);
boxplot(errs, 'Labels', string(N_dict_list));
set(gca, 'YScale', 'log');
xlabel('Dictionary size N_{dict}');
ylabel('Relative Frobenius reconstruction error (%, log scale)');
title(sprintf('Reconstruction Accuracy vs N_{dict}  (%d reps/point, %s)', R_reps, METHOD));
grid on;

%% -------------------- VIOLIN PLOT (self-contained, no toolbox needed) --------------------
figure('Position', [100, 100, 900, 500]);
hold on;

log_errs = log10(errs);                 % work in log-space since error spans decades
n_groups = numel(N_dict_list);
width_scale = 0.38;                     % max half-width of each violin

group_colors = lines(n_groups);

for k = 1:n_groups
    data_k = log_errs(:, k);
    data_k = data_k(~isnan(data_k));
    if numel(data_k) < 2, continue; end

    [f, xi] = ksdensity(data_k);        % f = density, xi = log10(error) support
    f = f / max(f) * width_scale;       % normalize violin half-width

    % mirrored fill: k - f(left side) ... k + f(right side)
    fill([k - f, fliplr(k + f)], [xi, fliplr(xi)], group_colors(k,:), ...
         'FaceAlpha', 0.5, 'EdgeColor', group_colors(k,:), 'LineWidth', 1.2);

    % overlay individual repetition points (jittered) and the median
    jitter = (rand(size(data_k)) - 0.5) * 0.12;
    plot(k + jitter, data_k, 'k.', 'MarkerSize', 8);
    plot([k-width_scale*0.5, k+width_scale*0.5], [median(data_k) median(data_k)], ...
         'k-', 'LineWidth', 2);
end

hold off;
set(gca, 'XTick', 1:n_groups, 'XTickLabel', string(N_dict_list));
xlabel('Dictionary size N_{dict}');

% Show y-axis in original % units even though violins are drawn in log10-space
yt = get(gca, 'YTick');
set(gca, 'YTickLabel', arrayfun(@(v) sprintf('%.3g', 10^v), yt, 'UniformOutput', false));
ylabel('Relative Frobenius reconstruction error (%)');
title(sprintf('Reconstruction Accuracy vs N_{dict}  (%d reps/point, %s) -- Violin', R_reps, METHOD));
grid on; xlim([0.5, n_groups+0.5]);

%% -------------------- Optional: local function for 'reseed' method --------------------

function [Xa_all, Ya_all] = generate_lorenz_pool(cfg)
    % cfg fields: dim, SCALE, SIGMA, BETA, RHO, M, delta_t, T_snap, T_burn
    % dim=3 reduces exactly to the original classic Lorenz case
    % (SCALE=1, matching the very first script in this conversation).
    dim = cfg.dim; SCALE = cfg.SCALE; SIGMA = cfg.SIGMA; BETA = cfg.BETA; RHO = cfg.RHO;
    M = cfg.M; delta_t = cfg.delta_t; T_snap = cfg.T_snap; T_burn = cfg.T_burn;

    asq  = 4/BETA - 1;
    didx = 1:20;
    dvec = ((2*didx-1).^2 + asq) / (1+asq);

    ODEFUN = build_generalized_lorenz_odefun(dim, SIGMA, BETA, RHO, SCALE, dvec);

    options = odeset('RelTol',1e-13,'AbsTol',1e-14);   % tighter tol, matches generalized-model runs

    X0_raw = 20*(rand(dim,M)-0.5);
    if dim >= 3
        X0_raw(3,:) = X0_raw(3,:) + 25;   % keep ICs near the attractor's "z" offset
    end

    Xa_ens = cell(M,1); Ya_ens = cell(M,1);
    for m = 1:M
        [~,Yb] = ode45(ODEFUN, [0 T_burn], X0_raw(:,m), options);
        x0 = Yb(end,:)';
        t_grid = 0:delta_t:(T_snap*delta_t);
        [~,Ytraj] = ode45(ODEFUN, t_grid, x0, options);
        Ytraj = Ytraj';
        Xa_ens{m} = Ytraj(:,1:end-1);
        Ya_ens{m} = Ytraj(:,2:end);
    end
    Xa_all = cell2mat(Xa_ens');
    Ya_all = cell2mat(Ya_ens');
end

function ODEFUN = build_generalized_lorenz_odefun(dim, SIGMA, BETA, RHO, SCALE, dvec)
    % Curbelo-style generalized Lorenz truncations (energy-conserving
    % quadratic extensions of the classic 3D Lorenz model).
    switch dim
        case 3
            ODEFUN = @(t,y) [SIGMA*(y(2)-y(1));
                y(1).*(RHO-y(3)*SCALE)-y(2);
                y(1).*y(2)*SCALE-BETA*y(3)];

        case 5
            ODEFUN = @(t,y) [SIGMA*(y(2)-y(1));
                y(1).*(RHO-y(3)*SCALE)-y(2);
                y(1).*y(2)*SCALE-BETA*y(3)-y(1).*y(4)*SCALE;
                y(1).*y(3)*SCALE-2*y(1).*y(5)*SCALE-dvec(2)*y(4);
                2*y(1).*y(4)*SCALE-4*BETA*y(5)];

        case 6
            ODEFUN = @(t,y) [SIGMA*(y(2)-y(1));
                -y(1).*y(3)*SCALE+y(4).*y(3)*SCALE-2*y(4).*y(6)*SCALE+RHO*y(1)-y(2);
                y(1).*y(2)*SCALE-y(1).*y(5)*SCALE-y(4).*y(2)*SCALE-BETA*y(3);
                -dvec(2)*SIGMA*y(4)+SIGMA/dvec(2)*y(5);
                y(1).*y(3)*SCALE-2*y(1).*y(6)*SCALE+RHO*y(4)-dvec(2)*y(5);
                2*y(1).*y(5)*SCALE+2*y(4).*y(2)*SCALE-4*BETA*y(6)];

        case 8
            ODEFUN = @(t,y) [SIGMA*(y(2)-y(1));
                -y(1).*y(3)*SCALE+y(4).*y(3)*SCALE-2*y(4).*y(6)*SCALE+RHO*y(1)-y(2);
                y(1).*y(2)*SCALE-BETA*y(3)-y(1).*y(5)*SCALE-y(4).*y(2)*SCALE-y(4).*y(7)*SCALE;
                -dvec(2)*SIGMA*y(4)+SIGMA/dvec(2)*y(5);
                y(1).*y(3)*SCALE-2*y(1).*y(6)*SCALE-dvec(2)*y(5)+RHO*y(4)-3*y(4).*y(8)*SCALE;
                2*y(1).*y(5)*SCALE-4*BETA*y(6)+2*y(4).*y(2)*SCALE-2*y(1).*y(7)*SCALE;
                2*y(1).*y(6)*SCALE-3*y(1).*y(8)*SCALE+y(4).*y(3)*SCALE-dvec(3)*y(7);
                3*y(1).*y(7)*SCALE+3*y(4).*y(5)*SCALE-9*BETA*y(8)];

        case 9
            ODEFUN = @(t,y) [SIGMA*(y(2)-y(1));
                -y(1).*y(3)*SCALE+RHO*y(1)-y(2)+y(4).*y(3)*SCALE-2*y(4).*y(6)*SCALE+2*y(7).*y(6)*SCALE-3*y(7).*y(9)*SCALE;
                y(1).*y(2)*SCALE-BETA*y(3)-y(1).*y(5)*SCALE-y(4).*y(2)*SCALE-y(4).*y(8)*SCALE-y(7).*y(5)*SCALE;
                -dvec(2)*SIGMA*y(4)+SIGMA/dvec(2)*y(5);
                y(1).*y(3)*SCALE-2*y(1).*y(6)*SCALE-dvec(2)*y(5)+RHO*y(4)-3*y(4).*y(9)*SCALE+y(7).*y(3)*SCALE;
                2*y(1).*y(5)*SCALE-4*BETA*y(6)+2*y(4).*y(2)*SCALE-2*y(1).*y(8)*SCALE-2*y(7).*y(2)*SCALE;
                -dvec(3)*SIGMA*y(7)+SIGMA/dvec(3)*y(8);
                2*y(1).*y(6)*SCALE-3*y(1).*y(9)*SCALE+y(4).*y(3)*SCALE-dvec(3)*y(8)+RHO*y(7);
                3*y(1).*y(8)*SCALE+3*y(4).*y(5)*SCALE-9*BETA*y(9)+3*y(7).*y(2)*SCALE];

        case 11
            ODEFUN = @(t,y) [SIGMA*(y(2)-y(1));
                -y(1).*y(3)*SCALE+RHO*y(1)-y(2)+y(4).*y(3)*SCALE-2*y(4).*y(6)*SCALE+2*y(7).*y(6)*SCALE-3*y(7).*y(9)*SCALE;
                y(1).*y(2)*SCALE-BETA*y(3)-y(1).*y(5)*SCALE-y(4).*y(2)*SCALE-y(4).*y(8)*SCALE-y(7).*y(5)*SCALE-y(7).*y(10)*SCALE;
                -dvec(2)*SIGMA*y(4)+SIGMA/dvec(2)*y(5);
                y(1).*y(3)*SCALE-2*y(1).*y(6)*SCALE-dvec(2)*y(5)+RHO*y(4)-3*y(4).*y(9)*SCALE+y(7).*y(3)*SCALE-4*y(7).*y(11)*SCALE;
                2*y(1).*y(5)*SCALE-4*BETA*y(6)+2*y(4).*y(2)*SCALE-2*y(1).*y(8)*SCALE-2*y(7).*y(2)*SCALE-2*y(4).*y(10)*SCALE;
                -dvec(3)*SIGMA*y(7)+SIGMA/dvec(3)*y(8);
                2*y(1).*y(6)*SCALE-3*y(1).*y(9)*SCALE+y(4).*y(3)*SCALE-dvec(3)*y(8)+RHO*y(7)-4*y(4).*y(11)*SCALE;
                3*y(1).*y(8)*SCALE+3*y(4).*y(5)*SCALE-9*BETA*y(9)+3*y(7).*y(2)*SCALE-3*y(1).*y(10)*SCALE;
                -4*y(1).*y(11)*SCALE+2*y(4).*y(6)*SCALE+3*y(1).*y(9)*SCALE+y(7).*y(3)*SCALE-dvec(4)*y(10);
                4*y(1).*y(10)*SCALE+4*y(4).*y(8)*SCALE+4*y(7).*y(5)*SCALE-16*BETA*y(11)];

        otherwise
            error('Unsupported dim=%d. Supported values: 3, 5, 6, 8, 9, 11.', dim);
    end
end