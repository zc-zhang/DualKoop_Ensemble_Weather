%% -------------------- BOXPLOT: Reconstruction Accuracy vs N_dict --------------------
% Assumes Xa_all, Ya_all already exist in the workspace (pooled Lorenz snapshots),
% and that kernel_ResDMD is on the path (your local version returning
% [G, K, L, PX, PY, PSI_x, PSI_y, PSI_y2, G1, A1, kernel_f]).

clc;
N_dict_list = [10 20 50 100 150 200 300];
R_reps      = 15;              % repetitions per N_dict value
frac_sub    = 0.85;            % fraction of pooled snapshots used per repetition

METHOD = 'bootstrap';          % 'bootstrap' (fast) or 'reseed' (regenerate Lorenz ensemble)

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
                [Xa_use, Ya_use] = generate_lorenz_pool();
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


%% -------------------- PLOT --------------------
figure('Position', [100, 100, 850, 480]);
boxplot(errs, 'Labels', string(N_dict_list));
set(gca, 'YScale', 'log');
xlabel('Dictionary size N_{dict}');
ylabel('Relative Frobenius reconstruction error (%, log scale)');
title(sprintf('Reconstruction Accuracy vs N_{dict}  (%d reps/point, %s)', R_reps, METHOD));
grid on;

%% -------------------- Optional: local function for 'reseed' method --------------------
function [Xa_all, Ya_all] = generate_lorenz_pool()
    SIGMA = 10; RHO = 28; BETA = 8/3;
    ODEFUN = @(t,y) [SIGMA*(y(2)-y(1));
                     y(1).*(RHO-y(3))-y(2);
                     y(1).*y(2)-BETA*y(3)];
    M = 50; 
    delta_t = 0.05; 
    T_snap = 40;
    T_burn = 5.0;
    options = odeset('RelTol',1e-10,'AbsTol',1e-11);
    X0_raw = 20*(rand(3,M)-0.5) + [0;0;25];
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


%% Additional testing by Residual using matrix L
