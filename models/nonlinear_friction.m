function [F_sx, F_fn, W_F_fn] = nonlinear_friction(dq_sx, r_smooth)
    % nonlinear_friction: Nonlinear regularized friction formulation
    import casadi.*
    n = length(dq_sx);
    p_fric = SX.sym('p_fric', 4*n);
    
    F_sx = SX.zeros(n,1);
    for i = 1:n
        fc_plus = p_fric(4*(i-1)+1);
        fc_minus = p_fric(4*(i-1)+2);
        fv_plus = p_fric(4*(i-1)+3);
        fv_minus = p_fric(4*(i-1)+4);
        
        v = dq_sx(i);
        t_sym = tanh(r_smooth * v);
        t_asym = (tanh(r_smooth * v)^2) / (abs(tanh(r_smooth * v)) + 1e-6);
        
        F_sym = 0.5*(fc_plus + fc_minus)*t_sym + 0.5*(fv_plus + fv_minus)*v;
        F_asym = 0.5*(fc_plus - fc_minus)*t_asym + 0.5*(fv_plus - fv_minus)*v*t_sym;
        
        F_sx(i) = F_sym + F_asym;
    end
    
    F_fn = Function('F_fn', {dq_sx, p_fric}, {F_sx});
    W_F_sx = jacobian(F_sx, p_fric);
    W_F_fn = Function('W_F_fn', {dq_sx}, {W_F_sx});
end
