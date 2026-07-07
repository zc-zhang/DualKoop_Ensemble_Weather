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