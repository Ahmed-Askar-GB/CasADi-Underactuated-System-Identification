function results = benchmark_solvers(W_D, tau_vec, filter_data, r_smooth)
    % benchmark_solvers: Computes Standard Unconstrained LS and Aguilar-Sierra Filtered LS
    
    if nargin < 4 || isempty(r_smooth)
        r_smooth = 100.0;
    end
    
    results = struct();
    
    %% 1. Standard Unconstrained Least Squares (Direct Regressor)
    fprintf('Running Standard Unconstrained Least Squares...\n');
    t_start_ls = tic;
    
    % Ridge regularization parameter for numerical stability if ill-conditioned
    lambda_ridge = 1e-8;
    Theta_ls = (W_D' * W_D + lambda_ridge * eye(size(W_D, 2))) \ (W_D' * tau_vec);
    time_ls = toc(t_start_ls);
    
    res_ls = W_D * Theta_ls - tau_vec;
    rmse_ls = sqrt(mean(res_ls.^2));
    sigma2_ls = mean(res_ls.^2);
    
    results.ls = struct();
    results.ls.Theta = Theta_ls;
    results.ls.residual = res_ls;
    results.ls.rmse = rmse_ls;
    results.ls.sigma2 = sigma2_ls;
    results.ls.time = time_ls;
    results.ls.status = 'converged';
    
    %% 2. Aguilar-Sierra Filtered Regressor LS (Baseline Method)
    if nargin >= 3 && ~isempty(filter_data)
        fprintf('Running Aguilar-Sierra Filtered Regressor Method...\n');
        t_start_ag = tic;
        
        l = filter_data.l;
        T = filter_data.T;
        qf1 = filter_data.qf1;
        qf2 = filter_data.qf2;
        dq1 = filter_data.dq1;
        dq2 = filter_data.dq2;
        u1 = filter_data.u1;
        u2 = filter_data.u2;
        
        r = r_smooth;
        
        % Omega_a matrix
        Oa11 = dq1;
        Oa12 = (dq1 .* sin(qf2)) .* sin(qf2);
        Oa13 = dq2 .* cos(qf2);
        Oa14 = zeros(l, 1); Oa15 = zeros(l, 1); Oa16 = zeros(l, 1);
        Oa17 = zeros(l, 1); Oa18 = zeros(l, 1); Oa19 = zeros(l, 1);

        Oa21 = zeros(l, 1); Oa22 = zeros(l, 1);
        Oa23 = dq1 .* cos(qf2);
        Oa24 = dq2;
        Oa25 = zeros(l, 1); Oa26 = zeros(l, 1); Oa27 = zeros(l, 1);
        Oa28 = zeros(l, 1); Oa29 = zeros(l, 1);

        % Omega_b matrix
        Ob11 = zeros(l, 1); Ob12 = zeros(l, 1); Ob13 = zeros(l, 1);
        Ob14 = zeros(l, 1); Ob15 = zeros(l, 1);
        Ob16 = dq1;
        Ob17 = zeros(l, 1);
        Ob18 = tanh(r * dq1);
        Ob19 = zeros(l, 1);

        Ob21 = zeros(l, 1);
        Ob22 = -0.5 * (dq1 .* dq1) .* sin(2 * qf2);
        Ob23 = (dq2 .* dq1) .* sin(qf2);
        Ob24 = zeros(l, 1);
        Ob25 = -sin(qf2);
        Ob26 = zeros(l, 1);
        Ob27 = dq2;
        Ob28 = zeros(l, 1);
        Ob29 = tanh(r * dq2);

        % Filter design
        Lambda = 30;
        B1 = [Lambda, -Lambda];
        A1 = [1, -exp(-Lambda * T)];
        B2 = [0, 1 - exp(-Lambda * T)];
        A2 = [1, -exp(-Lambda * T)];

        % Filter columns
        Oaf = zeros(l, 9, 2);
        Obf = zeros(l, 9, 2);
        
        Oa_cell = { {Oa11, Oa12, Oa13, Oa14, Oa15, Oa16, Oa17, Oa18, Oa19}, ...
                    {Oa21, Oa22, Oa23, Oa24, Oa25, Oa26, Oa27, Oa28, Oa29} };
        Ob_cell = { {Ob11, Ob12, Ob13, Ob14, Ob15, Ob16, Ob17, Ob18, Ob19}, ...
                    {Ob21, Ob22, Ob23, Ob24, Ob25, Ob26, Ob27, Ob28, Ob29} };
                
        for row = 1:2
            for col = 1:9
                Oaf(:, col, row) = filter(B1, A1, Oa_cell{row}{col});
                Obf(:, col, row) = filter(B2, A2, Ob_cell{row}{col});
            end
        end
        
        Of = Oaf + Obf;
        
        uf1 = filter(B2, A2, u1);
        uf2 = filter(B2, A2, u2);
        
        % Construct stacked filtered matrix
        Off = [Of(:,:,1); Of(:,:,2)];
        TauF = [uf1; uf2];
        
        Theta_ag = (Off' * Off) \ (Off' * TauF);
        time_ag = toc(t_start_ag);
        
        res_ag = Off * Theta_ag - TauF;
        rmse_ag = sqrt(mean(res_ag.^2));
        sigma2_ag = mean(res_ag.^2);
        
        results.aguilar = struct();
        results.aguilar.Theta = Theta_ag;
        results.aguilar.residual = res_ag;
        results.aguilar.rmse = rmse_ag;
        results.aguilar.sigma2 = sigma2_ag;
        results.aguilar.time = time_ag;
        results.aguilar.status = 'converged';
        results.aguilar.Off = Off;
        results.aguilar.TauF = TauF;
    end
end
