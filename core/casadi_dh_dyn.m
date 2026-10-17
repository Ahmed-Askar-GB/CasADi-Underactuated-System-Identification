function [M, C, G, N, Ro] = casadi_dh_dyn(DHpram, P_Ci, I, m, g, q, dq)
    % casadi_dh_dyn: CasADi SX-based port of DH_dyn
    import casadi.*
    n = length(q);
    P_Ci = P_Ci.'; % Transpose for consistency with original DH_dyn
    
    % Initialize transformation matrices
    Ro.L(1).A = Ai(DHpram(1,:));
    Ro.L(1).T0i = Ro.L(1).A;
    
    for i = 2:n
        Ro.L(i).A = Ai(DHpram(i,:));
        Ro.L(i).T0i = Ro.L(i-1).T0i * Ro.L(i).A;
    end
    
    % Initialize z-axis vectors
    z = SX.zeros(3, n);
    z(:,1) = [0;0;1];
    
    for i = 2:n
        z(:,i) = Ro.L(i-1).T0i(1:3,3);
    end
    
    % Compute CoM positions
    Ro.c(1).p = Ro.L(1).A(1:3,1:3) * P_Ci(:,1);
    for i = 2:n
        Ro.c(i).p = SX.zeros(3, i);
        Ro.c(i).p(:,1) = Ro.L(i-1).T0i(1:3,4) + Ro.L(i).T0i(1:3,1:3) * P_Ci(:,i);
        for j = 2:i
            Ro.c(i).p(:,j) = Ro.c(i).p(:,1) - Ro.L(j-1).T0i(1:3,4);
        end
    end
    
    % Compute Jacobians
    for i = 1:n
        Ro.L(i).jv = SX.zeros(3, n);
        Ro.L(i).jw = SX.zeros(3, n);
        for j = 1:n
            if j <= i
                is_revolute = (full(evalf(DHpram(j,5))) == 1);
                if is_revolute % Revolute joint
                    Ro.L(i).jv(:,j) = cross(z(:,j), Ro.c(i).p(:,j));
                    Ro.L(i).jw(:,j) = z(:,j);
                else % Prismatic joint
                    Ro.L(i).jv(:,j) = z(:,j);
                    Ro.L(i).jw(:,j) = SX.zeros(3,1);
                end
            end
        end
    end
    
    % Gravity vector
    G = SX.zeros(n,1);
    for i = 1:n
        for j = 1:n
            G(i) = G(i) - m(j) * g.' * Ro.L(j).jv(:,i);
        end
    end
    
    % Inertia matrix
    M = SX.zeros(n,n);
    for i = 1:n
        if iscell(I)
            I_i = I{i};
        else
            I_i = I(:,:,i);
        end
        I0i = Ro.L(i).T0i(1:3,1:3) * I_i * Ro.L(i).T0i(1:3,1:3).';
        M = M + m(i)*(Ro.L(i).jv.' * Ro.L(i).jv) + Ro.L(i).jw.' * I0i * Ro.L(i).jw;
    end
    
    % Coriolis matrix (Christoffel symbols via exact Jacobians)
    C = SX.zeros(n,n);
    M_vec = M(:);
    dM_dq = jacobian(M_vec, q); % n^2 x n
    
    for i = 1:n
        for j = 1:n
            for k = 1:n
                dMij_dqk = dM_dq((j-1)*n + i, k);
                dMik_dqj = dM_dq((k-1)*n + i, j);
                dMkj_dqi = dM_dq((j-1)*n + k, i);
                C(i,j) = C(i,j) + 0.5 * (dMij_dqk + dMik_dqj - dMkj_dqi) * dq(k);
            end
        end
    end
    
    % Passivity matrix N = dM/dt - 2C
    Mdot_vec = jtimes(M_vec, q, dq);
    Mdot = reshape(Mdot_vec, n, n);
    N = Mdot - 2*C;
end

function A = Ai(Dh)
    % A_i transformation matrix
    ai = Dh(1); alphai = Dh(2); si = Dh(3); di = Dh(4);
    A = [cos(si), -sin(si)*cos(alphai),  sin(si)*sin(alphai), ai*cos(si); 
         sin(si),  cos(si)*cos(alphai), -cos(si)*sin(alphai), ai*sin(si); 
         0,        sin(alphai),          cos(alphai),         di; 
         0,        0,                    0,                   1];
end
