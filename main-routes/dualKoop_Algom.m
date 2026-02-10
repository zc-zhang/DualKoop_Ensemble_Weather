% Another lok at the dual Koopman opertor in paper 
% see the Algorithm 1 in paper 
% J. Mohet, A. Mauroy, and J. Winkin, A dual Koopman approach to observer design for nonlinear systems, 
% arXiv;2503.08345v1.    % https://arxiv.org/pdf/2503.08345

function [G,A,K_star,lambda,W,C,PSI_x,PSI_y,PSI_y2,kernel_f, UU] ...
    = dualKoop_Algom(Xa,Ya,varargin)
% =========================================================================
% Dual Koopman operator via kernelized Residual DMD
% (Algorithm 1 compatible, RKHS adjoint Koopman)
%
% OUTPUTS (Algorithm 1 notation emphasized):
%   G        : Gram matrix <k_xi, k_xj>
%   A        : Cross Gram matrix <k_yi, k_xj>
%   K_star   : Dual Koopman operator (whitened)
%   lambda   : Koopman eigenvalues
%   W        : Eigenvectors in orthonormal RKHS basis
%   C        : Dual Koopman modes c_j (Algorithm 1 coefficients)
%   PSI_*    : Evaluations of eigenfunctions on data
% =========================================================================

% ---------------------------
% Parse inputs
% ---------------------------
p = inputParser;
addParameter(p,'N',size(Xa,2))
addParameter(p,'type',"Gaussian")
addParameter(p,'cut_off',0)
addParameter(p,'Xb',[])
addParameter(p,'Yb',[])
addParameter(p,'Y2',[])
addParameter(p,'UU',[],@isnumeric);
parse(p,varargin{:})

% ---------------------------
% Kernel definitions
% ---------------------------
d = mean(vecnorm(Xa-mean(Xa,2)));

if isnumeric(p.Results.type)
    % Polynomial
    kernel_f = @(x,y) (y'*x/d^2 + 1).^p.Results.type;
    
elseif p.Results.type=="Linear"
    kernel_f = @(x,y) y'*x;
    
elseif p.Results.type=="Gaussian"
    % Corrected Gaussian kernel
    kernel_f = @(x,y) exp(-(vecnorm(x).^2' + vecnorm(y).^2 - 2*real(y'*x)) / d^2);
    
elseif p.Results.type=="Laplacian"
    kernel_f = @(x,y) exp(-sqrt(vecnorm(x).^2' + vecnorm(y).^2 - 2*real(y'*x)) / d);
    
elseif p.Results.type=="Lorentzian"
    kernel_f = @(x,y) (1 + (vecnorm(x).^2' + vecnorm(y).^2 - 2*real(y'*x)) / d^2).^(-1);
    
elseif p.Results.type=="Cauchy"
    kernel_f = @(x,y) 1 ./ (1 + (vecnorm(x).^2' + vecnorm(y).^2 - 2*real(y'*x)) / d^2);
    
else
    error("Unknown kernel type")
end
% ---------------------------
% Algorithm 1: Gram matrices
% ---------------------------
G1 = kernel_f(Xa,Xa);   G1 = (G1+G1')/2;
A1 = kernel_f(Ya,Xa)';  
L1 = kernel_f(Ya,Ya);   L1 = (L1+L1')/2;
% ---------------------------
% RKHS whitening (G^{-1/2})
% ---------------------------
[U,D] = eig(G1 + norm(G1)*p.Results.cut_off*eye(size(G1)));

[d_sorted,idx] = sort(real(diag(D)),'descend');
U = U(:,idx);
D = diag(d_sorted);

N = min(p.Results.N, sum(d_sorted > 0));
U = U(:,1:N);
D = D(1:N,1:N);

if isempty(p.Results.UU)
    UU = U * diag(1 ./ sqrt(diag(D)));   % G^{-1/2}
else
    UU = p.Results.UU;
end



% Algorithm 1 objects (in orthonormal basis)
G = UU'*G1*UU;           % ≈ I
A = UU'*A1*UU;
K_star = A;              % dual Koopman operator

% ---------------------------
% Eigendecomposition
% ---------------------------
[W,Dk] = eig(K_star);
lambda = diag(Dk);

% ---------------------------
% Dual Koopman modes (IMPORTANT)
% c_j = G^{-1/2} w_jwhat n -
% ---------------------------
C = UU * W;

% ---------------------------
% Eigenfunction evaluation
% phi_j(x) = sum_i c_ij k(x,x_i)
% ---------------------------
PSI_x  = G1 * C;   % Evaluate on Xa itself
PSI_y  = A1' * C;  % Evaluate on Ya itself
PSI_y2 = [];

if ~isempty(p.Results.Xb)
    PSI_x = kernel_f(p.Results.Xb,Xa)' * C;
end
if ~isempty(p.Results.Yb)
    PSI_y = kernel_f(p.Results.Yb,Xa)' * C;
end
if ~isempty(p.Results.Y2)
    PSI_y2 = kernel_f(p.Results.Y2,Xa)' * C;
end

end
