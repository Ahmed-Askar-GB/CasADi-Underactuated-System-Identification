function constraints = extract_constraints(M_sx, N_sx, p_sx, q_sx, dq_sx, q_pts)
    % extract_constraints: Automated constraints for positive-definiteness & skew-symmetry
    % Extracted via Algorithmic Differentiation over the Computational Graph.
    %
    % Addresses Reviewer 1 & Reviewer 3:
    % 1. Evaluates row gradients across independent configurations q_pts and diverse
    %    velocity samples, isolating the true linearly independent row space via SVD.
    % 2. Automatically computes all leading principal minors Delta_1, ..., Delta_n
    %    for general n-DOF mechanisms (including n=2 RIFP and n=3 DIPC).
    import casadi.*
    n = size(M_sx, 1);
    np = length(p_sx);
    
    if nargin < 6 || isempty(q_pts)
        % Default diverse sampling spanning reachable workspace envelope
        q_pts = [zeros(1, n);
                 (pi/4)*eye(n);
                 -(pi/4)*eye(n);
                 (pi/2)*ones(1, n);
                 -(pi/3)*ones(1, n)];
    end
    Nc = size(q_pts, 1);

    % --- C1: Symmetry Equality Constraints (D_ij(q) = D_ji(q)) ---
    A_sym_raw = [];
    for i = 1:n
        for j = i+1:n
            diff_sx = M_sx(i,j) - M_sx(j,i);
            row = jacobian(diff_sx, p_sx);
            row_fn = Function('row_fn', {q_sx}, {row});
            for k = 1:Nc
                q_k = q_pts(k,:)';
                A_sym_raw = [A_sym_raw; full(row_fn(q_k))];
            end
        end
    end
    % Extract linearly independent rows via SVD
    if ~isempty(A_sym_raw)
        [~, S, V] = svd(A_sym_raw, 'econ');
        diag_s = diag(S);
        tol = max(size(A_sym_raw)) * eps(max(diag_s));
        r_sym = sum(diag_s > max(tol, 1e-9));
        if r_sym > 0
            constraints.A_sym = V(:, 1:r_sym)';
        else
            constraints.A_sym = zeros(0, np);
        end
    else
        constraints.A_sym = zeros(0, np);
    end

    % --- C3: Skew-Symmetry Passivity Constraints (N_ij + N_ji = 0) ---
    A_skew_raw = [];
    % Generate non-zero velocity test points
    dq_test = [ones(1, n); -ones(1, n); 2*eye(n); -2*eye(n)];
    N_vel = size(dq_test, 1);
    for i = 1:n
        for j = i+1:n
            skew_expr = N_sx(i,j) + N_sx(j,i);
            row = jacobian(skew_expr, p_sx);
            row_fn = Function('skew_row_fn', {q_sx, dq_sx}, {row});
            for k = 1:min(Nc, 10)
                q_k = q_pts(k,:)';
                for v = 1:N_vel
                    dq_v = dq_test(v,:)';
                    A_skew_raw = [A_skew_raw; full(row_fn(q_k, dq_v))];
                end
            end
        end
    end
    if ~isempty(A_skew_raw)
        [~, S, V] = svd(A_skew_raw, 'econ');
        diag_s = diag(S);
        tol = max(size(A_skew_raw)) * eps(max(diag_s));
        r_skew = sum(diag_s > max(tol, 1e-9));
        if r_skew > 0
            constraints.A_skew = V(:, 1:r_skew)';
        else
            constraints.A_skew = zeros(0, np);
        end
    else
        constraints.A_skew = zeros(0, np);
    end

    % --- C2: Positive Definiteness: Leading Principal Minors ---
    % Evaluates Delta_1, ..., Delta_n at each configuration point
    constraints.minors_fns = cell(Nc, n);
    constraints.D11_fns = cell(Nc, 1);
    constraints.detM_fns = cell(Nc, 1);

    for m = 1:n
        sub_M = M_sx(1:m, 1:m);
        det_sub = det(sub_M);
        for k = 1:Nc
            q_k = q_pts(k,:)';
            det_sub_k = substitute(det_sub, q_sx, SX(q_k));
            fn_k = Function(sprintf('minor_%d_%d', m, k), {p_sx}, {det_sub_k});
            constraints.minors_fns{k, m} = fn_k;
            if m == 1
                constraints.D11_fns{k} = fn_k;
            end
            if m == n
                constraints.detM_fns{k} = fn_k;
            end
        end
    end
end
