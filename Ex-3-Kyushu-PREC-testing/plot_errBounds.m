%% plot_error_bound.m

%% --- 1. Setup Variables and Matrices ---
KEs_sorted  = KEs(sort_idx);
KEFs_sorted = KEFs(:, sort_idx);
%KMs_sorted  = KMs(:, sort_idx);

r_values = [50,100,150,200,250,300]; 
num_r = length(r_values);
pinv_tol = 1e-2;

fprintf('size(Ya_all)  = [%d, %d]\n', size(Ya_all));
fprintf('size(Ya_mean) = [%d, %d]\n', size(Ya_mean));
total_snapshots = size(Ya_all, 1);      % IN-SAMPLE: snapshots in training ensemble

%% ---- Storage ----
log_residual_in  = zeros(total_snapshots, num_r);   % IN-SAMPLE  (per training snapshot)
log_residual_out = zeros(num_r, 1);                  % OUT-OF-SAMPLE (single target: Ya_mean)
err_bound_curve  = zeros(num_r, 1);                  % theoretical bound, same for both panels

total_eigen_mass = sum(abs(KEs_sorted));

for r_idx = 1:num_r
    r = r_values(r_idx);
    KEs_r        = KEs_sorted(1:r);
    KEFs_r       = KEFs_sorted(:, 1:r);
    KEFs_mean_r  = KEFs_mean_sorted(:, 1:r);

    KMs_r = Xa_all * pinv(KEFs_r.', pinv_tol);   % fit once per r, on training data ONLY

    % ---- IN-SAMPLE: Xa_all -> Ya_all ----
    Ya_pred_r = real(KMs_r * diag(KEs_r) * KEFs_r.');
    true_norms  = sqrt(sum(Ya_all.^2, 2));
    error_norms = sqrt(sum((Ya_all - Ya_pred_r).^2, 2));
    res_in = error_norms ./ true_norms;
    res_in(res_in < 1e-12) = 1e-12;
    log_residual_in(:, r_idx) = log10(res_in);

    % ---- OUT-OF-SAMPLE: Xa_mean -> Ya_mean (same KMs_r, same KEs_r) ----
    Ya_pred_mean_r = real(KMs_r * diag(KEs_r) * KEFs_mean_r.');
    res_out = norm(Ya_mean - Ya_pred_mean_r, 'fro') / norm(Ya_mean, 'fro');
    if res_out < 1e-12, res_out = 1e-12; end
    log_residual_out(r_idx) = log10(res_out);

    % ---- Theoretical bound (same formula, shared reference for both) ----
    unmodeled_eigen_mass = sum(abs(KEs_sorted(r+1:end)));
    bound_val = unmodeled_eigen_mass / total_eigen_mass;
    if bound_val < 1e-12, bound_val = 1e-12; end
    err_bound_curve(r_idx) = log10(bound_val);
end

%% ---- Plot: two panels, in-sample (with error bars) | out-of-sample (single curve) ----
res_in_mean = mean(log_residual_in, 1);
res_in_std  = std(log_residual_in, 0, 1);

fig = figure('Position', [100, 100, 1100, 450]);

subplot(1,2,1);
errorbar(r_values, res_in_mean, res_in_std, '-s', ...
    'Color',[0.2 0.4 0.8], 'LineWidth',1.8, 'MarkerSize',7, ...
    'MarkerFaceColor',[0.2 0.4 0.8], 'CapSize',6); hold on;
plot(r_values, err_bound_curve, '--o', 'Color',[0.8 0.2 0.2], 'LineWidth',1.8, 'MarkerSize',6);
grid on;
xlabel('Retained Modes (r)','FontWeight','bold');
ylabel('log_{10}(Relative Error)','FontWeight','bold');
title('In-Sample (Training Fit)','FontWeight','bold');
legend({'Empirical (mean \pm std over snapshots)','Theoretical Bound'}, 'Location','best');

subplot(1,2,2);
plot(r_values, log_residual_out, '-o', 'Color',[0.85 0.30 0.25], 'LineWidth',1.8, ...
    'MarkerSize',7, 'MarkerFaceColor',[0.85 0.30 0.25]); hold on;
plot(r_values, err_bound_curve, '--o', 'Color',[0.8 0.2 0.2], 'LineWidth',1.8, 'MarkerSize',6);
grid on;
xlabel('Retained Modes (r)','FontWeight','bold');
ylabel('log_{10}(Relative Error)','FontWeight','bold');
title('Out-of-Sample (Ensemble Mean Forecast)','FontWeight','bold');
legend({'Empirical','Theoretical Bound'}, 'Location','best');

sgtitle('In-Sample vs. Out-of-Sample: Empirical Error vs. Theoretical Bound', 'FontWeight','bold');