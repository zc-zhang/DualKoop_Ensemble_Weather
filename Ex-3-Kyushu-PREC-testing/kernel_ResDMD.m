% Dual Koopman operator can be computed by means of Residual DMD
%  Github by M. Colbrook
% https://github.com/MColbrook/Residual-Dynamic-Mode-Decomposition/tree/main/main_routines
%      


%function [G,K_star,L,PX,PY,PSI_x,PSI_y,PSI_y2] = kernel_ResDMD(Xa,Ya,varargin)
%function [G,K_star,L,PX,PY,PSI_x,PSI_y,PSI_y2,G1,A1,kernel_f,KEs,KMs,KEFs] = kernel_ResDMD(Xa,Ya,varargin)
function [G,K,L,PX,PY,PSI_x,PSI_y,PSI_y2,G1,A1,UU, kernel_f] = kernel_ResDMD(Xa,Ya,varargin)
% This code applies kernelized ResDMD.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% Xa and Ya: data matrices used in kernel_EDMD to form dictionary (columns
% correspond to instances of the state variable)

% OPTIONAL LABELLED INPUTS
% N: size of computed dictionary, default is number of data points for kernel EDMD
% type: kernel used, default is normalised Gaussian, "Laplacian" is for
% nomralised Laplacian, and numeric value (e.g., 20) is for polynomial
% kernel
% cut_off: stability parameter for SVD, default is 0
% Xb, Yb: additional data matrices used in ResDMD for test data
% Y2: additional data matrix for stochastic version

% OUTPUTS
% G K L matrices for kernelResDMD
% PSI matrices for ResDMD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect the optional inputs
p = inputParser;

addParameter(p,'N',size(Xa,2),@(x) x==floor(x))
addParameter(p,'type',"Gaussian");
addParameter(p,'cut_off',0,@(x) x>=0)
addParameter(p,'Xb',[],@isnumeric)
addParameter(p,'Yb',[],@isnumeric)
addParameter(p,'Y2',[],@isnumeric)
addParameter(p,'c',[],@isnumeric)       % NEW: polynomial scaling factor
addParameter(p,'alpha',[],@isnumeric)   % NEW: polynomial power

p.CaseSensitive = false;
parse(p,varargin{:})

% Apply kernel EDMD
if isnumeric(p.Results.type)
    d = mean(vecnorm(Xa));
    kernel_f = @(x,y) (y'*x/d^2+1).^(p.Results.type);
elseif p.Results.type=="Linear"
    kernel_f = @(x,y) y'*x;
elseif p.Results.type=="Polynomial"
    if isempty(p.Results.c)
        c = mean(vecnorm(Xa));   % default scaling, same convention as other kernels
    else
        c = p.Results.c;
    end
    if isempty(p.Results.alpha)
        error('Polynomial kernel requires alpha to be specified, e.g. ''alpha'', 3')
    end
    kernel_f = @(x,y) (y'*x/c^2 + 1).^p.Results.alpha;
elseif p.Results.type=="Laplacian"
    d = mean(vecnorm(Xa-mean(Xa,2)));
    if isa(Xa,'single') % safeguard against square root (but a little bit slower)
        kernel_f = @(x,y) exp(-pdist2(y',x')/d);
    else
        kernel_f = @(x,y) exp(-sqrt(-2*real(y'*x)+dot(x,x)+dot(y,y)')/d);
    end
elseif p.Results.type=="Gaussian"
    d = mean(vecnorm(Xa-mean(Xa,2)));
    kernel_f = @(x,y) exp(-(-2*real(y'*x)+dot(x,x)+dot(y,y)')/d^2);
elseif p.Results.type=="Lorentzian"
    d = mean(vecnorm(Xa-mean(Xa,2)));
    kernel_f = @(x,y) (1+(-2*real(y'*x)+dot(x,x)+dot(y,y)')/d^2).^(-1);
end

% (i) Kernel formulation -> infinite-dimensional RKHS
G1 = kernel_f(Xa,Xa); 
G1 = (G1+G1')/2;
A1 = kernel_f(Ya,Xa)'; % dual without transpose 
%A1_naive = kernel_f(Ya,Xa)'; % naive 
L1 = kernel_f(Ya,Ya);  L1 = (L1+L1')/2;

% Post processing
% (ii) SVD + cutoff = residual-norm control
[U,D0] = eig(G1+norm(G1)*p.Results.cut_off*eye(size(G1)));

[~,I] = sort(diag(D0),'descend');
U = U(:,I); D0 = D0(I,I);
N = min(p.Results.N,length(find(diag(D0)>0)));
U = U(:,1:N); D0 = D0(1:N,1:N);
UU = U*sqrt(diag(1./diag(D0)));

 G = UU'*G1*UU; % approximtae I_N, % i.e., G = eye(N);
% (iii) orthonormalization in RKHS norm
%==============================================
%---- K = G1^{-1/2}AG1^{-1/2}  
% -> This is not the EDMD opertaor G^{-1}A   -> symmetric, norm-consistent, residual-optimal in RKHS.
% -> Compuet the dual Koopman operator K
% --> In kernel EDMD. Kernelized ResEDMD 
% ---> basis function are kernel sections k_{x_i}, 
%  evolution g_t+1 =K*g_t    (coefficients live in the dual space)
%==================================== =========
K = UU'*A1*UU; %dual % KeDMD of KeDMD
%K=K';

% K_naive =UU'A1_naive*UU  % K_naive = K_satr'


% stochastic ResDMD for operator L 
L = UU'*L1*UU;

%   here K_star is essnetiall meets 
%  <K*g,h> =<g,Kh>, \forall g,h \in span{k_{x_i}}

PX = G1'*UU;
PY = A1*UU;

if ~isempty(p.Results.Xb) % test data case
    PSI_x = kernel_f(p.Results.Xb,Xa)'*UU;
else
    PSI_x =[];
end

if ~isempty(p.Results.Yb) % test data case
    PSI_y = kernel_f(p.Results.Yb,Xa)'*UU;
else
    PSI_y =[];
end

if ~isempty(p.Results.Y2) % stochastic case
    PSI_y2 = kernel_f(p.Results.Y2,Xa)'*UU;
else
    PSI_y2 =[];
end


end