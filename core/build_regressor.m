function [Y_sx, Y_min_fn, Theta_def] = build_regressor(M, C, G, q, dq, ddq, p_sx)
    % build_regressor: Symbolic-to-numeric CasADi regressor generator
    import casadi.*
    n = length(q);
    np = length(p_sx);
    
    % EOM residual (linear in parameters)
    tau_sx = M * ddq + C * dq + G;
    
    % Regressor matrix via CasADi Jacobian
    Y_sx = jacobian(tau_sx, p_sx);
    
    % QR Base Reduction
    Y_fn = Function('Y_fn', {q, dq, ddq, p_sx}, {Y_sx});
    
    N_qr = max(4*np, 200);
    W_stack = zeros(N_qr*n, np);
    
    % Populate W_stack with random numeric evaluations
    for k = 1:N_qr
        q_k = rand(n,1)*2*pi - pi;
        dq_k = randn(n,1)*5;
        ddq_k = randn(n,1)*20;
        W_stack((k-1)*n+1 : k*n, :) = full(Y_fn(q_k, dq_k, ddq_k, zeros(np,1)));
    end
    
    % QR Decomposition
    [~, R, E] = qr(W_stack, 'vector');
    diag_R = abs(diag(R));
    tol = max(size(W_stack)) * eps(max(diag_R)) * 1000;
    r = sum(diag_R > tol);
    
    idx_base = E(1:r);
    idx_dep = E(r+1:end);
    
    R1 = R(1:r, 1:r);
    R2 = R(1:r, r+1:end);
    beta = R1 \ R2;
    beta(abs(beta) < 1e-10) = 0;
    
    Y_min_sx = Y_sx(:, idx_base);
    Y_min_fn = Function('Y_min', {q, dq, ddq}, {Y_min_sx});
    
    Theta_def.idx_base = idx_base;
    Theta_def.idx_dep = idx_dep;
    Theta_def.beta = beta;
    Theta_def.r = r;
end
