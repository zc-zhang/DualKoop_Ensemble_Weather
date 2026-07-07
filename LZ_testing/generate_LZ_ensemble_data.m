function [Xa,Ya,w,xQ] = generate_LZ_ensemble_data(dim,N,delta_t,SCALE)
% GENERATE_LZ_ENSEMBLE_DATA  Ensemble (quadrature-node) snapshot data for
% the generalized Lorenz (LZ) hierarchy, suitable for kernel (Res)DMD.
%
% This reproduces the data-generation logic of your original script, but
% packaged as a function so you can loop over dim = 3,5,6,8,9,11 and feed
% the result straight into kernel_ResDMD / weighted_kernel_ResDMD.
%
% Key idea: instead of one long trajectory, we sample M initial conditions
% x^(m) on a Gauss-Hermite sparse grid (nodes = xQ, weights = w, w>=0,
% sum(w)=1), and propagate each ONE step forward by delta_t to get
% y^(m) = F_{delta_t}(x^(m)). {x^(m), w^(m), y^(m)} is precisely the
% "quadrature data" ResDMD wants (Colbrook & Townsend 2024, Sec 2.2):
% it approximates integrals against a reference measure (here, the
% Gaussian reference measure from Hermite quadrature) without ever
% needing ergodicity / a single long orbit.
%
% INPUTS
%   dim     : 3,5,6,8,9, or 11 (which generalized Lorenz model)
%   N       : sparse-grid resolution parameter (your script used N=16)
%   delta_t : forward propagation time (your script used 0.05)
%   SCALE   : scaling parameter in the ODE (your script used 1/10)
%
% OUTPUTS
%   Xa : dim x M matrix of initial conditions (quadrature nodes)
%   Ya : dim x M matrix of one-step-forward images, Ya(:,m)=F(Xa(:,m))
%   w  : M x 1 vector of quadrature weights, w>=0, sum(w)=1
%   xQ : M x dim (same info as Xa', kept for convenience/plotting)

if nargin<1, dim=3; end
if nargin<2, N=5; end  % less than 16
if nargin<3, delta_t=0.05; end
if nargin<4, SCALE=1/10; end

SIGMA=10; BETA=8/3; RHO=40;
asq=4/BETA-1;
dvec=1:20;
dvec=((2*dvec-1).^2+asq)/(1+asq);

switch dim
    case 3
        ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));
                        y(1).*(RHO-y(3)*SCALE)-y(2);
                        y(1).*y(2)*SCALE-BETA*y(3)];
    case 5
        ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));
                        y(1).*(RHO-y(3)*SCALE)-y(2);
                        y(1).*y(2)*SCALE-BETA*y(3)-y(1).*y(4)*SCALE;
                        y(1).*y(3)*SCALE-2*y(1).*y(5)*SCALE-dvec(2)*y(4);
                        2*y(1).*y(4)*SCALE-4*BETA*y(5)];
    case 6
        ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));
                        -y(1).*y(3)*SCALE+y(4).*y(3)*SCALE-2*y(4).*y(6)*SCALE+RHO*y(1)-y(2);
                        y(1).*y(2)*SCALE-y(1).*y(5)*SCALE-y(4).*y(2)*SCALE-BETA*y(3);
                        -dvec(2)*SIGMA*y(4)+SIGMA/dvec(2)*y(5);
                        y(1).*y(3)*SCALE-2*y(1).*y(6)*SCALE+RHO*y(4)-dvec(2)*y(5);
                        2*y(1).*y(5)*SCALE+2*y(4).*y(2)*SCALE-4*BETA*y(6)];
    case 8
        ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));
                        -y(1).*y(3)*SCALE+y(4).*y(3)*SCALE-2*y(4).*y(6)*SCALE+RHO*y(1)-y(2);
                        y(1).*y(2)*SCALE-BETA*y(3)-y(1).*y(5)*SCALE-y(4).*y(2)*SCALE-y(4).*y(7)*SCALE;
                        -dvec(2)*SIGMA*y(4)+SIGMA/dvec(2)*y(5);
                        y(1).*y(3)*SCALE-2*y(1).*y(6)*SCALE-dvec(2)*y(5)+RHO*y(4)-3*y(4).*y(8)*SCALE;
                        2*y(1).*y(5)*SCALE-4*BETA*y(6)+2*y(4).*y(2)*SCALE-2*y(1).*y(7)*SCALE;
                        2*y(1).*y(6)*SCALE-3*y(1).*y(8)*SCALE+y(4).*y(3)*SCALE-dvec(3)*y(7);
                        3*y(1).*y(7)*SCALE+3*y(4).*y(5)*SCALE-9*BETA*y(8)];
    case 9
        ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));
                        -y(1).*y(3)*SCALE+RHO*y(1)-y(2)+y(4).*y(3)*SCALE-2*y(4).*y(6)*SCALE+2*y(7).*y(6)*SCALE-3*y(7).*y(9)*SCALE;
                        y(1).*y(2)*SCALE-BETA*y(3)-y(1).*y(5)*SCALE-y(4).*y(2)*SCALE-y(4).*y(8)*SCALE-y(7).*y(5)*SCALE;
                        -dvec(2)*SIGMA*y(4)+SIGMA/dvec(2)*y(5);
                        y(1).*y(3)*SCALE-2*y(1).*y(6)*SCALE-dvec(2)*y(5)+RHO*y(4)-3*y(4).*y(9)*SCALE+y(7).*y(3)*SCALE;
                        2*y(1).*y(5)*SCALE-4*BETA*y(6)+2*y(4).*y(2)*SCALE-2*y(1).*y(8)*SCALE-2*y(7).*y(2)*SCALE;
                        -dvec(3)*SIGMA*y(7)+SIGMA/dvec(3)*y(8);
                        2*y(1).*y(6)*SCALE-3*y(1).*y(9)*SCALE+y(4).*y(3)*SCALE-dvec(3)*y(8)+RHO*y(7);
                        3*y(1).*y(8)*SCALE+3*y(4).*y(5)*SCALE-9*BETA*y(9)+3*y(7).*y(2)*SCALE];
    case 11
        ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));
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
        error('dim=%d not implemented (use 3,5,6,8,9,11)',dim);
end

options = odeset('RelTol',1e-13,'AbsTol',1e-14);


%% Sparse-grid (quadrature) nodes = ensemble of initial conditions
warning('off','SparseGKit:uint16')
[lev2knots,idxset]=define_functions_for_rule('SM',dim);
knots=@(n) hermpts(n);
[S,~]= create_sparse_grid(dim,log2(2*N)+4,knots,lev2knots,idxset); % <-- BROKEN LINE
Sr=reduce_sparse_grid(S);

xQ=transpose(Sr.knots);                    % M x dim nodes
w =transpose(Sr.weights)/sqrt(pi^dim);      % M x 1 weights, sum(w)~1

M=size(xQ,1);
Ya_=zeros(M,dim);

 pf = parfor_progress(M);
 pfcleanup = onCleanup(@() delete(pf));
parfor j=1:M
    Y0=transpose(xQ(j,:));
    [~,Y]=ode45(ODEFUN,[0.000001 delta_t 2*delta_t],Y0,options);
    Ya_(j,:)=Y(2,:);
    parfor_progress(pf);
end

Xa = xQ';   % dim x M  (kernel_ResDMD convention: columns = data points)
Ya = Ya_';  % dim x M
w  = w(:);

% --- SAVE THE DATA HERE ---
% Creates a unique filename based on the dimension, e.g., 'LZ_data_dim5.mat'
filename = sprintf('LZ_ensemble_data_dim%d.mat', dim);
save(filename, 'Xa', 'Ya', 'w', 'xQ');
fprintf('Successfully generated and saved %d points to %s\n', M, filename);
end