%% Plot_ACC

% %% =========================================================================
% %  ACC COMPUTATION FOR RESHAPED ENSEMBLE DATA (20860 x 7000)  -- corrected
% %  =========================================================================
% d      = 20860; % Spatial points (p)
% M      = 100;   % Number of ensemble members / cross-validated cases
% N_snap = 70;    % Number of time snapshots
% 
% O_tensor = reshape(Ya_all,  [d, N_snap, M]);
% O_tensor = permute(O_tensor, [1, 3, 2]);   % [20860 x 100 x 70]
% F_tensor = reshape(Ya_pred, [d, N_snap, M]);
% F_tensor = permute(F_tensor, [1, 3, 2]);   % [20860 x 100 x 70]
% 
% ACC_matrix = zeros(M, N_snap);
% for t = 1:N_snap
%     O_t = O_tensor(:, :, t);   % [20860 x 100]
%     F_t = F_tensor(:, :, t);   % [20860 x 100]
% 
%     % --- SHARED climatology/reference (truth-based), same for O and F ---
%     C_t    = mean(O_t, 2);           % [20860 x 1]
%     O_anom = O_t - C_t;              % [20860 x 100]
%     F_anom = F_t - C_t;              % <-- fixed: subtract the SAME reference
% 
%     for i = 1:M
%         num = sum(O_anom(:, i) .* F_anom(:, i));
%         den = sqrt(sum(O_anom(:, i).^2)) * sqrt(sum(F_anom(:, i).^2));
%         ACC_matrix(i, t) = (den > 0) * (num / max(den, eps));
%     end
% end
% 
% ACC_mean = mean(ACC_matrix, 1);
% ACC_std  = std(ACC_matrix, 0, 1);
% upper_bound = min(1, ACC_mean + ACC_std);
% lower_bound = max(-1, ACC_mean - ACC_std);
% t_axis = 1:N_snap;
% 
% figure('Color', 'w', 'Position', [150, 150, 750, 450]);
% hold on;
% fill([t_axis, fliplr(t_axis)], [upper_bound, fliplr(lower_bound)], ...
%      [0.75, 0.85, 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
%      'DisplayName', 'Ensemble Spread (\pm1 Std)');
% plot(t_axis, ACC_matrix', 'Color', [0.7, 0.7, 0.7, 0.2], ...
%      'LineWidth', 0.6, 'HandleVisibility', 'off');
% plot(t_axis, ACC_mean, 'b-', 'LineWidth', 2.5, ...
%      'DisplayName', 'Dual KMD Ensemble Mean ACC');
% yline(0.6, 'r--', 'Skill Cutoff (0.6)', 'LineWidth', 1.5, ...
%       'DisplayName', '0.6 Skill Threshold');
% xlabel('Time Snapshot Step (t)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Anomaly Correlation Coefficient (ACC)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Kyushu PREC Dual KMD ACC Skill Across 100 Ensemble Members', ...
%       'FontSize', 13, 'FontWeight', 'bold');
% ylim([-1, 1.05]);   % don't clip decay/negative values while diagnosing
% xlim([1, N_snap]);
% grid on; box on;
% legend('Location', 'southwest', 'FontSize', 10);


%=======================
%% ============================================================
%% Per-member spatial ACC(t) — gives spread across M=100 members
%% ============================================================
% Ya_raw_3D, Ya_pred_3D are [20860 x 70 x 100]  (already reshaped earlier)

nT = size(Ya_raw_3D, 2);   % 70
M  = size(Ya_raw_3D, 3);   % 100

EnACC_members = nan(M, nT);   % [100 x 70]

% Reference climatology: event-mean of the RAW ensemble-mean field
% (same reference used for all members, for consistency)
clim_ref = mean(Ya_raw_ens_mean, 2);   % [20860 x 1]

for k = 1:M
    Anom_raw_k  = Ya_raw_3D(:,:,k)  - clim_ref;   % [20860 x 70]
    Anom_pred_k = Ya_pred_3D(:,:,k) - clim_ref;   % [20860 x 70]

    for t = 1:nT
        f = Anom_pred_k(:, t);
        o = Anom_raw_k(:, t);
        EnACC_members(k, t) = sum(f.*o) / ( sqrt(sum(f.^2))*sqrt(sum(o.^2)) );
    end
end


%% Ensemble-mean-based EnACC(t)  (already derived earlier — repeated here)
Anom_raw_EM  = Ya_raw_ens_mean  - clim_ref;
Anom_pred_EM = Ya_pred_ens_mean - clim_ref;

EnACC_em = nan(1, nT);
for t = 1:nT
    f = Anom_pred_EM(:, t);
    o = Anom_raw_EM(:, t);
    EnACC_em(t) = sum(f.*o) / ( sqrt(sum(f.^2))*sqrt(sum(o.^2)) );
end



%% ============================================================
%% FIGURE: EnACC spread (shading) vs Ensemble-Mean EnACC (line)
%% ============================================================
prc_lo = prctile(EnACC_members, 5,  1);   % [1 x 70]  5th percentile across members
prc_hi = prctile(EnACC_members, 95, 1);   % [1 x 70]  95th percentile
EnACC_avg_members = mean(EnACC_members, 1, 'omitnan');  % member-averaged ACC (EnACC_avg)

figure('Color','w','Position',[100 100 1000 550]);
hold on;

% Shaded spread region (5th-95th percentile across 100 members)
fill([timeVec, fliplr(timeVec)], [prc_lo, fliplr(prc_hi)], ...
     [0.7 0.85 1], 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
     'DisplayName', 'Member spread (5–95%)');

% Member-averaged ACC (mean of the 100 individual ACC values)
plot(timeVec, EnACC_avg_members, 'b--', 'LineWidth', 1.5, ...
     'DisplayName', 'EnACC_{avg} (mean of members)');

% Ensemble-mean-based ACC (ACC on the averaged field itself)
plot(timeVec, EnACC_em, 'r-', 'LineWidth', 2.5, ...
     'DisplayName', 'EnACC_{em} (ensemble-mean field)');

hold off; grid on; ylim([-1 1]);
xlabel('Date / Time (UTC)'); ylabel('EnACC');
title('Spatial EnACC: Spread (shading) vs Ensemble-Mean Skill');
legend('Location','southwest');

ax = gca;
ax.XTick = tickTimes;
set(ax, 'XTickLabel', tickLabels, 'XTickLabelRotation', 25, 'FontSize', 11);


%% Kumamoto: member spread (shading) + ensemble-mean lines
raw_Kumamoto_members  = squeeze(Ya_raw_3D(Kumamoto_Locatidx, :, :));   % [70 x 100]
pred_Kumamoto_members = squeeze(Ya_pred_3D(Kumamoto_Locatidx, :, :));  % [70 x 100]

raw_lo  = prctile(raw_Kumamoto_members,  5, 2)';  raw_hi  = prctile(raw_Kumamoto_members,  95, 2)';
pred_lo = prctile(pred_Kumamoto_members, 5, 2)';  pred_hi = prctile(pred_Kumamoto_members, 95, 2)';

figure('Color','w'); hold on;
fill([timeVec,fliplr(timeVec)], [raw_lo,fliplr(raw_hi)],  [0.6 1 0.6], 'FaceAlpha',0.3,'EdgeColor','none');
fill([timeVec,fliplr(timeVec)], [pred_lo,fliplr(pred_hi)],[1 0.7 1], 'FaceAlpha',0.3,'EdgeColor','none');
plot(timeVec, data_raw_EM_Kumamoto,  'g-',  'LineWidth',2,'DisplayName','Raw Ensemble Mean');
plot(timeVec, data_pred_EM_Kumamoto, 'm--', 'LineWidth',2,'DisplayName','Reconstructed Ensemble Mean');
hold off; grid on; legend('Location','northwest');
title('Kumamoto: Member Spread (shading) + Ensemble-Mean Comparison');



%% ============================================================
%% Kumamoto Point: EnACC(t) — Member Spread (shading) + 
%%                 Ensemble-Mean ACC (line)
%% y-axis = ACC (unitless), using windowed correlation
%% ============================================================

W  = 9;     % moving-window half-width in hours (25h window); tune as needed
nT = 70;
M  = 100;

% --- Reference/"climatology" at Kumamoto: event time-mean of raw EM ---
clim_ref_p = mean(data_raw_EM_Kumamoto);   % scalar

% --- Extract all-member trajectories at Kumamoto ---
Kumamoto_raw_allmembers  = squeeze(Ya_raw_3D(Kumamoto_Locatidx, :, :));   % [70 x 100]
Kumamoto_pred_allmembers = squeeze(Ya_pred_3D(Kumamoto_Locatidx, :, :));  % [70 x 100]

% --- Per-member windowed ACC(t) ---
ACC_Kumamoto_members = nan(M, nT);   % [100 x 70]
for k = 1:M
    o_full = Kumamoto_raw_allmembers(:,k)'  - clim_ref_p;   % [1x70]
    f_full = Kumamoto_pred_allmembers(:,k)' - clim_ref_p;   % [1x70]
    for t = 1:nT
        idx = max(1,t-W):min(nT,t+W);
        o = o_full(idx);  f = f_full(idx);
        ACC_Kumamoto_members(k,t) = sum(f.*o) / ( sqrt(sum(f.^2))*sqrt(sum(o.^2)) );
    end
end

% --- Ensemble-mean-based windowed ACC(t) (EnACC_em) ---
o_EM = data_raw_EM_Kumamoto  - clim_ref_p;
f_EM = data_pred_EM_Kumamoto - clim_ref_p;
EnACC_Kumamoto_em = nan(1,nT);
for t = 1:nT
    idx = max(1,t-W):min(nT,t+W);
    o = o_EM(idx); f = f_EM(idx);
    EnACC_Kumamoto_em(t) = sum(f.*o) / ( sqrt(sum(f.^2))*sqrt(sum(o.^2)) );
end

% --- Member-averaged ACC (mean of the 100 individual ACC curves) ---
EnACC_Kumamoto_avg = mean(ACC_Kumamoto_members, 1, 'omitnan');

% --- Spread bounds across members (5th-95th percentile) ---
ACC_lo = prctile(ACC_Kumamoto_members, 5,  1);   % [1 x 70]
ACC_hi = prctile(ACC_Kumamoto_members, 95, 1);   % [1 x 70]

%% ---- Plot: ACC as the sole y-axis quantity ----
figure('Color','w','Name','Kumamoto EnACC','Position',[100 100 1000 550]);
hold on;

% Shaded spread of per-member ACC
fill([timeVec, fliplr(timeVec)], [ACC_lo, fliplr(ACC_hi)], ...
     [0.7 0.85 1], 'FaceAlpha', 0.45, 'EdgeColor', 'none', ...
     'DisplayName', 'Member ACC spread (5–95%)');

% Member-averaged ACC
plot(timeVec, EnACC_Kumamoto_avg, 'b--', 'LineWidth', 1.5, ...
     'DisplayName', 'EnACC_{avg} (mean of members)');

% Ensemble-mean-based ACC
plot(timeVec, EnACC_Kumamoto_em, 'r-', 'LineWidth', 2.5, ...
     'DisplayName', 'EnACC_{em} (ensemble-mean field)');

% Reference line at 0 (no skill)
yline(0, 'k:', 'LineWidth', 1);

hold off; grid on; ylim([-1 1]);
xlabel('Date / Time (UTC)', 'FontSize', 12);
ylabel('ACC', 'FontSize', 12);
title(sprintf('Point EnACC at Kumamoto (%d-hr moving window)', 2*W+1), ...
      'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 10);

ax = gca;
ax.XTick = tickTimes;
set(ax, 'XTickLabel', tickLabels, 'XTickLabelRotation', 25, 'FontSize', 11);