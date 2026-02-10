% plot _dual _Koopman _using 'kernel_ResDMD.m'

% plot_Kernel_ResDMD.m

clear; clc; close all;

%% ============================================================
% 1. Generate Lorenz data
%% ============================================================

sigma = 10; beta = 8/3; rho = 50;
dt = 0.01; T = 25;
tspan = 0:dt:T;

x0 = [1;1;1];
[~,X] = ode45(@(t,x) lorenz63(t,x,sigma,beta,rho), tspan, x0);
X = X';                         % 3 × Nt

Xa = X(:,1:end-1);
Ya = X(:,2:end);

%% ============================================================
% 2. Kernel ResDMD (dual Koopman in RKHS)
%% ============================================================

[G,K,L,PX,PY,PSI_x] = kernel_ResDMD( ...
    Xa, Ya, ...
    'type','Gaussian', ...
    'N',300, ...
    'Xb', Xa ...
);

% option 2:
%[G,K,L,PX,PY,PSI_x] = kernel_ResDMD(Xa,Ya,'type','Gaussian','N',300);


[W,D] = eig(K);
lambda = diag(D);

%% ============================================================
% 3. Dual Koopman eigenvalues
%% ============================================================

figure;
plot(real(lambda), imag(lambda), 'ko','MarkerFaceColor','k');
xlabel('Re(\lambda)');
ylabel('Im(\lambda)');
%title('Dual Koopman eigenvalues (kernel ResDMD)');
title('Dual Koopman eigenvalues');
grid on; axis equal;

%% ============================================================
% 4. Dual Koopman modes (kernel coefficients)
%% ============================================================

figure;
stem(abs(W(:,2)),'filled');
xlabel('RKHS mode index');
ylabel('|Coefficient|');
title('Dual Koopman mode (eigenvector 2)');
grid on;

%% ============================================================
% 5. Dual Koopman eigenfunction (Lorenz butterfly)
%% ============================================================

% Use PSI_x from kernel_ResDMD
% PSI_x : Nt × N  (feature evaluation)

eig_id = 2;
% PSI_x = kernel_f(Xtest,Xa)' * UU;
% phi   = PSI_x * W(:,eig_id);

phi = real(PSI_x * W(:,eig_id));   % Nt × 1

figure;
scatter3(Xa(1,:), Xa(2,:), Xa(3,:), ...
         25, phi, 'filled');
axis equal
view(3)
colorbar
title(['Dual Koopman eigenfunction \phi_',num2str(eig_id)]);
xlabel('x'); ylabel('y'); zlabel('z');

%% ============================================================
% Lorenz system
%% ============================================================
function dx = lorenz63(~,x,sigma,beta,rho)
dx = zeros(3,1);
dx(1) = sigma*(x(2) - x(1));
dx(2) = x(1)*(rho - x(3)) - x(2);
dx(3) = x(1)*x(2) - beta*x(3);
end
