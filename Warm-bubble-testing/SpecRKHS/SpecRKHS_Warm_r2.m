%% Koopman/PF analysis of Northern hemisphere sea level heights

%% load data sets
%clear
rng(0)


%% =====================================================================
Xa_all = EnVfull_Xa_36;
Ya_all = EnVfull_Ya_36;
p      = size(Xa_all, 1);
M      = 10;
T_snap = size(Xa_all, 2) / M;
grid_size = [97, 40];
nLAND = (1:p).';

T_train_split = 29;                   % First 29 columns for training
T_len =29;
steps =T_snap-T_len;

X_train_list  = cell(1, M);
Y_train_list  = cell(1, M);
X_test_list   = cell(1, M);
Y_test_list   = cell(1, M);
%Kx0_vals_cell = cell(1, M);

for m = 1:M
    col_start = (m - 1) * T_snap + 1;
    
    % 1. Training data (First 29 columns)
    X_train_list{m} = Xa_all(:, col_start : col_start + T_train_split - 1);
    Y_train_list{m} = Ya_all(:, col_start : col_start + T_train_split - 1);
    
    % 2. Testing data (Remaining columns from split point to end)
    test_start_idx  = col_start + T_train_split - 1;
    test_end_idx    = col_start + T_snap - 1;
    X_test_list{m}  = Xa_all(:, test_start_idx : test_end_idx);
    Y_test_list{m}  = Ya_all(:, test_start_idx : test_end_idx);
end

% Combine training sets across all 10 ensembles (Size: 3880 x 290)
X_train = horzcat(X_train_list{:}); 
Y_train = horzcat(Y_train_list{:});

% Combine test sets across all 10 ensembles
X_test  = horzcat(X_test_list{:});
Y_test  = horzcat(Y_test_list{:});
x=X_train;
y=Y_train;

%% 
% --- Optional: Define Ensembles cell array if you need it ---
Ensembles = cell(1, M);
for m = 1:M
    col_start = (m - 1) * T_snap + 1;
    col_end   = m * T_snap;
    Ensembles{m} = Xa_all(:, col_start:col_end);
end
%% compute matrices
ker=@(x,t) kernel_matern(x,t);
[G,A,R]=generate_matrices_kernelized(x,y,ker);

%% compute verified eigenvalues
num=290;
[Lambda_res,F_res,Lambda,F,res,res_verif,idx,W,W_res]=verified_eigenvalues(G,A,R,num);
length(idx)

%% compute Perron-Frobenius modes and plot
for idx=[2 4]
    L=Lambda_res(idx); 
    F=F_res(:,idx); 
    r=res_verif(idx);
    Phi=((G*F)\(x.')).';
    
    figure
    u = Phi;
    %u = real(u*exp(1i*mean(angle(u))));
    v = zeros(97*40,1)+NaN;
    v(nLAND) = u(:);
    v = reshape(v,[97,40]);
    %v=flip(v.');
    %imagesc(v,'AlphaData',~isnan(v))
    ax = nexttile;
    imagesc(data.y, data.z, abs(v), 'AlphaData', ~isnan(v))
    colormap(ax, brighten(redblueTecplot(21), -0.55));
    colorbar;
    %clim([mean(u(:)) - 2*std(u(:)), mean(u(:)) + 2*std(u(:))])
    set(gca, 'Color', [1,1,1]*0.6)
    axis xy; axis equal
    xlim([0 2e4]); ylim([0 2e4])
    set(gca, 'FontSize', 10);
    xlabel("y", 'FontSize', 10); ylabel("z", 'FontSize', 10);
   % title(sprintf('$\\textrm{Mode }%d: \\lambda_{%d} = %.2f %+.2f i$', j, j, real(L_j), imag(L_j)), ...
    %    'Interpreter', 'latex', 'FontSize', 10);
    hold on;
    plot(data.y(26), data.z(25), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow');
    hold off;
end

%% plot spurious and verified eigenvalues
figure
scatter(angle(Lambda),log(abs(Lambda)),200,res,'.','LineWidth',1);
hold on
scatter(angle(Lambda_res),log(abs(Lambda_res)),500,res_verif,'.','LineWidth',1);
box on
clim([0,0.1])
load('cmap.mat')
colormap(cmap2); colorbar
xlabel('$\mathrm{arg}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{log}(|\lambda|)$','interpreter','latex','fontsize',18)
title(['Residuals for sea level data',newline],'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18; axis([-pi pi -20*10^(-3) 10^(-3)])
for k=-5:1:5
    plot(k*pi/6*ones(22,1),-20*10^(-3):10^(-3):10^(-3),'--','Color','black')
end
xticks([-pi -5*pi/6 -4*pi/6 -3*pi/6 -2*pi/6 -pi/6 0 pi/6 2*pi/6 3*pi/6 4*pi/6 5*pi/6 pi])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
xticklabels({'$-\pi$','','$-2\pi/3$','','$-\pi/3$','','$0$','','$\pi/3$','','$2\pi/3$','','$\pi$'})
xtickangle(30)
exportgraphics(gcf,'sea_level_evals_angle.pdf','ContentType','vector','BackgroundColor','none')

%% compute psuedospectra
% pts=100;
% x_pts=linspace(-1.2,1.2,pts);    y_pts=linspace(-0.02,1.2,pts/2);
% z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);
% res_pspec=pseudospectra(G,A,R,z_pts);

%% plot pseudospectral contours
% res_pspec_rs=reshape(res_pspec,length(y_pts),length(x_pts));
% figure
% hold on
% box on
% v=(10.^(-10:0.2:0));
% contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(res_pspec_rs)),log10(v));
% contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(res_pspec_rs)),log10(v));
% cbh=colorbar;
% cbh.Ticks=log10(10.^(-2:1:0));
% cbh.TickLabels=10.^(-2:1:0);
% clim([-2,0]);
% reset(gcf)
% set(gca,'YDir','normal')
% colormap(inferno(100))
% axis equal;
% 
% plot(sin(0:0.01:2*pi), cos(0:0.01:2*pi),'--','color','white'); %plot unit circle
% 
% title('Pseudospectrum for sea level data','interpreter','latex','fontsize',18)
% xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
% ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
% 
% ax=gca; ax.FontSize=18; axis equal tight;   axis([x_pts(1),x_pts(end),-y_pts(end),y_pts(end)])
% exportgraphics(gcf,'sea_level_pspec.pdf','ContentType','vector','BackgroundColor','none')

%% =====================================================================
%  COMPUTE PREDICTIONS & ERRORS ACROSS ALL ENSEMBLES (Fully Corrected)
%  =====================================================================
num_train = size(X_train, 2);
% Pre-allocate error tracking matrices for all ensembles and lead times
% Size: [M ensembles x steps lead times]
err_kmd_all   = zeros(M, steps);
err_kedmd_all = zeros(M, steps);
err_dmd_all   = zeros(M, steps);

%% --- Pre-calculate Standard DMD globally (since it doesn't depend on individual x0 for the operator) ---
[U_dmd, S_dmd, ~] = svd(X_train, 'econ');
r = rank(S_dmd);
U_dmd = U_dmd(:, 1:r);
PXs = X_train' * U_dmd;
PYs = Y_train' * U_dmd;
K_dmd = PXs \ PYs;
[W_dmd, LAM_dmd] = eig(K_dmd);
lambda_dmd_vec = diag(LAM_dmd);
PXr = PXs * W_dmd; 

for m = 1:M
    data_m = Ensembles{m};
    
    % 1. Define initial state x0 and extract ground truth for the test window
    x0        = data_m(:, T_train_split);
    real_data = data_m(:, T_train_split + (1:steps)); % Ground truth future states [States x steps]
    
    %% --- A. SpecRKHS-Obs Predictions ---
    Kx0_vals = zeros(num_train, 1);
    for i = 1:num_train
        Kx0_vals(i) = ker(x0, X_train(:, i));
    end
    
    coefs_res = ((G * F_res) \ Kx0_vals); 
    time_factor_kmd = (Lambda_res(:) .^ (1:steps)); 
    temporal_dynamics_kmd = coefs_res(:) .* time_factor_kmd; % [num_res x steps]
    x_kmd = real( (X_train * F_res) * temporal_dynamics_kmd );        
    
    %% --- B. kEDMD Predictions ---
    G_start = zeros(1, num_train);
    for i = 1:num_train
        G_start(i) = ker(x0, X_train(:, i));
    end
    
    mode_full = (([G; G_start] * W) \ (Y_train.')).'; 
    psi0_full = G_start * W;                          
    
    time_factor_kedmd = (conj(Lambda(:)) .^ (1:steps)); % [num_train x steps]
    temporal_dynamics_kedmd = psi0_full(:) .* time_factor_kedmd; % [num_train x steps]
    x_kedmd = real( mode_full * temporal_dynamics_kedmd ); 
    
    %% --- C. Standard DMD Predictions (Corrected) ---
    % Project individual initial state x0 into the reduced DMD subspace
    b_dmd = PXr \ (U_dmd' * x0);
    
    % Temporal evolution in the reduced subspace: [r x steps]
    time_factor_dmd = (lambda_dmd_vec .^ (1:steps));
    subspace_dynamics = b_dmd .* time_factor_dmd; 
    
    % Map back to full physical space: [States x r] * [r x steps] -> [States x steps]
    x_dmd = real( (U_dmd * W_dmd) * subspace_dynamics );
    
    %% --- D. Calculate Relative Forecast Errors for this Ensemble ---
    err_kmd_all(m, :)   = sum(abs(x_kmd   - real_data).^2, 1) ./ sum(abs(real_data).^2, 1);
    err_kedmd_all(m, :) = sum(abs(x_kedmd - real_data).^2, 1) ./ sum(abs(real_data).^2, 1);
    err_dmd_all(m, :)   = sum(abs(x_dmd   - real_data).^2, 1) ./ sum(abs(real_data).^2, 1);
end
%=========================================
%  AVERAGE ERRORS ACROSS ALL ENSEMBLES & PLOT
%  =====================================================================
er1 = mean(err_dmd_all, 1);    % Mean DMD error
er2 = mean(err_kedmd_all, 1);  % Mean kEDMD error
er3 = mean(err_kmd_all, 1);    % Mean SpecRKHS-Obs error

figure;
plot(er1, 'linewidth', 2)
hold on
plot(er2, 'linewidth', 2)
plot(er3, 'linewidth', 2)
grid on

title('Relative forecast errors comparison (Ensemble Average)', 'fontsize', 18, 'interpreter', 'latex')
xlabel('Lead Time (Steps)', 'interpreter', 'latex', 'fontsize', 18)
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18)
legend({'DMD', 'kEDMD', 'SpecRKHS-Obs'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'best')
ax = gca; ax.FontSize = 18; box on;

exportgraphics(gcf, 'sea_level_error_delay.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

%% =====================================================================
%  WEATHER ENSEMBLE OBSERVABLE FORECAST (Ensemble Averaged)
%  =====================================================================
num_train = size(X_train, 2);

% Pre-allocate storage for tracking across M ensembles
weather_exact_all   = zeros(M, steps + 1);
weather_kmd_all     = zeros(M, steps + 1);
weather_kedmd_all   = zeros(M, steps + 1);

for m = 1:M
    data_m = Ensembles{m};
    t_train_len = size(X_train_list{m}, 2);
    
    % 1. Observable history from training portion (spatial mean across 3880 grid points)
    weather_old = mean(data_m(:, 1:t_train_len), 1).'; % [t_train_len x 1]
    
    % 2. SpecRKHS-Obs Mode computation for this ensemble
    mode_res = (G * F_res) \ weather_old;
    
    % Initial state and future real test window for this ensemble
    x0          = data_m(:, T_train_split);
    real_window = data_m(:, T_train_split + (0:steps)); % [3880 x (steps+1)]
    
    % Kernel vector for x0 against X_train
    Kx0_vals = zeros(num_train, 1);
    for i = 1:num_train
        Kx0_vals(i) = ker(x0, X_train(:, i));
    end
    coefs_res = ((G * F_res) \ Kx0_vals).';
    
    % --- SpecRKHS-Obs Forecast of Weather Observable ---
    weather_kmd_all(m, :) = real((coefs_res .* ((Lambda_res).^(0:steps)).') * (F_res.' * G * F_res * mode_res));
    
    % --- kEDMD Forecast of Weather Observable ---
    weather_old_full = mean(data_m(:, 1:(t_train_len + 1)), 1).'; 
    G_start = zeros(1, num_train);
    for i = 1:num_train
        G_start(i) = ker(x0, X_train(:, i));
    end
    psi0_full = G_start * W;
    mode_full = (([G; G_start]*W) \ weather_old_full).';
    weather_kedmd_all(m, :) = real(transpose(transpose(psi0_full) .* (conj(Lambda).^(0:steps))) * mode_full.')';
    
    % --- Exact Values (Spatial mean of the true weather field over time) ---
    weather_exact_all(m, :) = mean(real_window, 1); 
end

%% =====================================================================
%  AVERAGE ACROSS ALL ENSEMBLES & PLOT
%  =====================================================================
weather_exact   = mean(weather_exact_all, 1);
weather_kmd     = mean(weather_kmd_all, 1);
weather_kedmd   = mean(weather_kedmd_all, 1);

figure;
p1 = plot(0:1:steps, weather_kedmd, 'linewidth', 2, 'color', [0.8500 0.3250 0.0980]);
hold on;
p2 = plot(0:1:steps, weather_kmd,   'linewidth', 2, 'color', [0.9290 0.6940 0.1250]);
p3 = plot(0:1:steps, weather_exact, '--', 'linewidth', 2, 'color', [0 0.4470 0.7410]);

grid on;
title('Weather Ensemble Field Forecast Comparison', 'fontsize', 18, 'interpreter', 'latex');
xlabel('Lead Time (Steps)', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Spatial Mean Observable Value', 'interpreter', 'latex', 'fontsize', 18);
legend([p3 p1 p2], {'Exact', 'kEDMD', 'SpecRKHS-Obs'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'northeast');

ax = gca; ax.FontSize = 18; box on;
exportgraphics(gcf, 'weather_observable_prediction_mean.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%% =====================================================================
%  FIGURE 1: Weather Forecast Without kEDMD (Exact vs. SpecRKHS-Obs)
%  =====================================================================
% Plot comparing the exact spatial mean weather path against 
% the SpecRKHS-Obs prediction over the forecast lead steps.
figure;
hold on;
p2 = plot(0:1:steps, weather_kmd,   'linewidth', 2, 'color', [0.9290 0.6940 0.1250]);
p3 = plot(0:1:steps, weather_exact, 'linewidth', 2, 'color', [0 0.4470 0.7410]);
grid on;
title('Weather Ensemble Forecast Comparison', 'fontsize', 18, 'interpreter', 'latex');
xlabel('Lead Time (Steps)', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Spatial Mean Weather Value', 'interpreter', 'latex', 'fontsize', 18);
legend([p3 p2], {'Exact', 'SpecRKHS-Obs'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'southwest');
exportgraphics(gcf, 'weather_prediction_nokedmd.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');


%% =====================================================================
%  FIGURE 2: Relative Error of SpecRKHS-Obs vs Observable Exact Values
%  =====================================================================
% Semilog plot of the relative error between the exact observable 
% function values and the SpecRKHS-Obs prediction.
figure;
er_obs = (abs(weather_observable_exact - weather_kmd).' .^ 2) ./ (abs(weather_observable_exact) .^ 2);
semilogy(0:steps, er_obs, 'linewidth', 2, 'color', [0.9290 0.6940 0.1250]);
grid on;
title('Weather Observable Relative Forecast Error', 'fontsize', 18, 'interpreter', 'latex');
xlabel('Lead Time (Steps)', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18);
legend({'SpecRKHS-Obs'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'northeast');
exportgraphics(gcf, 'weather_error_observable.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');


%% =====================================================================
%  FIGURE 3: Relative Error of SpecRKHS-Obs vs Exact Mean Weather Values
%  =====================================================================
% Semilog plot tracking the relative forecasting error 
% between the true ensemble mean and the SpecRKHS-Obs projection.
figure;
er_mean = (abs(weather_exact - weather_kmd).^2) ./ (abs(weather_exact).^2);
semilogy(0:1:steps, er_mean, 'linewidth', 2, 'color', [0.9290 0.6940 0.1250]);
grid on;
title('Weather Ensemble Relative Forecast Error', 'fontsize', 18, 'interpreter', 'latex');
xlabel('Lead Time (Steps)', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18);
legend({'SpecRKHS-Obs'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'northeast');
exportgraphics(gcf, 'weather_error_mean.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');


%% =====================================================================
%  FIGURE 4: Three-Way Comparison (Exact, Observable, and SpecRKHS-Obs)
%  =====================================================================
% Multi-curve comparison showing the exact mean weather field, 
% the exact observable tracking, and the SpecRKHS-Obs prediction together.
figure;
hold on;
p3 = plot(0:1:steps, weather_exact, 'linewidth', 2, 'color', [0 0.4470 0.7410]);
p1 = plot(0:1:steps, weather_observable_exact, 'linewidth', 2, 'color', [0.4940 0.1840 0.5560]);
p2 = plot(0:1:steps, weather_kmd, 'linewidth', 2, 'color', [0.9290 0.6940 0.1250]);
grid on;
title('Weather Forecast Comparison (Observable vs Exact)', 'fontsize', 18, 'interpreter', 'latex');
xlabel('Lead Time (Steps)', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Weather Field Value', 'interpreter', 'latex', 'fontsize', 18);
legend([p3 p1 p2], {'Exact', 'Observable', 'SpecRKHS-Obs'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'southwest');
exportgraphics(gcf, 'weather_prediction_observable.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');


%% =====================================================================
%  FIGURE 5: Kernel Function $K_{x_0}$ Relative Forecast Error
%  =====================================================================
% Semilog plot showing the relative error behavior specifically 
% for the RKHS kernel evaluation $K_{x_0}$ over time.
figure;
er_kernel = (abs(Kx0_future_vals - Kx0_predict).^2) ./ (abs(Kx0_future_vals).^2);
semilogy(0:steps, er_kernel, 'linewidth', 2, 'color', [0.9290 0.6940 0.1250]);
grid on;
title('Relative Forecast Error for Kernel Function', 'fontsize', 18, 'interpreter', 'latex');
xlabel('Lead Time (Steps)', 'interpreter', 'latex', 'fontsize', 18);
ylabel('Relative Forecast Error', 'interpreter', 'latex', 'fontsize', 18);
legend({'SpecRKHS-Obs'}, 'interpreter', 'latex', 'fontsize', 16, 'location', 'northeast');
exportgraphics(gcf, 'weather_prediction_kernel_error.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%% matern kernel d=60330, n-d/2=2
function ker=kernel_matern(x,t)
    sigma=1/10000;
    r=vecnorm(x-t);
    ker=zeros(1,size(x,2));
    ker(r>0)=(sigma*r(r>0)).^2.*besselk(-2,sigma*r(r>0));
    ker(r==0)=2;
end