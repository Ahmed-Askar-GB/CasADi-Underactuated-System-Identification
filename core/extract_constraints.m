function constraints = extract_constraints(M_sx, N_sx, p_sx, q_sx, dq_sx, q_pts)
    % extract_constraints: Automated constraints for positive-definiteness & skew-symmetry
    import casadi.*
    n = size(M_sx,1);
    
    % C1: Symmetry constraints
    A_sym = [];
    for i = 1:n
        for j = i+1:n
            diff_sx = M_sx(i,j) - M_sx(j,i);
            row = jacobian(diff_sx, p_sx);
            row_fn = Function('row_fn', {q_sx}, {row});
            A_sym = [A_sym; full(row_fn(zeros(n,1)))];
        end
    end
    constraints.A_sym = A_sym;
    
    % C3: Skew-symmetry constraints
    A_skew = [];
    for i = 1:n
        for j = i+1:n
            skew_expr = N_sx(i,j) + N_sx(j,i);
            row = jacobian(skew_expr, p_sx);
            row_fn = Function('row_fn', {q_sx, dq_sx}, {row});
            % Evaluated at zero state to get constant linear mapping
            A_skew = [A_skew; full(row_fn(zeros(n,1), zeros(n,1)))];
        end
    end
    constraints.A_skew = A_skew;
    
    % C2: Positive Definiteness constraints (Pointwise evaluation)
    Nc = size(q_pts, 1);
    constraints.D11_fns = cell(Nc, 1);
    constraints.detM_fns = cell(Nc, 1);
    
    M11 = M_sx(1,1);
    M12 = M_sx(1,2);
    M22 = M_sx(2,2);
    detM = M11*M22 - M12^2;
    
    for k = 1:Nc
        q_k = q_pts(k,:)';
        M11_k = substitute(M11, q_sx, SX(q_k));
        detM_k = substitute(detM, q_sx, SX(q_k));
        
        constraints.D11_fns{k} = Function(sprintf('D11_%d',k), {p_sx}, {M11_k});
        constraints.detM_fns{k} = Function(sprintf('detM_%d',k), {p_sx}, {detM_k});
    end
end
