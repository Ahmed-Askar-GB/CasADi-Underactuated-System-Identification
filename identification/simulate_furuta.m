function [t_sim, q_sim, dq_sim, sim_metrics] = simulate_furuta(th, t, u, q_meas, dq_meas, r_smooth)
    % simulate_furuta: Forward dynamic simulation of the Furuta pendulum using RK4
    % Inputs:
    %   th: 9x1 parameter vector
    %   t: time vector (Nx1)
    %   u: motor torque vector (Nx1)
    %   q_meas: measured joint angles (Nx2, optional)
    %   dq_meas: measured joint velocities (Nx2, optional)
    %   r_smooth: friction smoothing factor (default 100.0)
    
    if nargin < 6 || isempty(r_smooth)
        r_smooth = 100.0;
    end
    
    N = length(t);
    T = mean(diff(t));
    
    if nargin >= 4 && ~isempty(q_meas)
        x0 = [q_meas(1,1); q_meas(1,2); dq_meas(1,1); dq_meas(1,2)];
    else
        x0 = [0; pi; 0; 0];
    end
    
    x = zeros(4, N);
    x(:, 1) = x0;
    
    for k = 1:N-1
        u_k = u(k);
        
        % RK4 integration steps
        k1 = furuta_deriv(x(:, k), u_k, th, r_smooth);
        k2 = furuta_deriv(x(:, k) + 0.5 * T * k1, u_k, th, r_smooth);
        k3 = furuta_deriv(x(:, k) + 0.5 * T * k2, u_k, th, r_smooth);
        k4 = furuta_deriv(x(:, k) + T * k3, u_k, th, r_smooth);
        
        x(:, k+1) = x(:, k) + (T / 6) * (k1 + 2*k2 + 2*k3 + k4);
    end
    
    t_sim = t;
    q_sim = [x(1, :)', x(2, :)'];
    dq_sim = [x(3, :)', x(4, :)'];
    
    sim_metrics = struct();
    if nargin >= 4 && ~isempty(q_meas)
        sim_metrics.rmse_q1 = sqrt(mean((q_sim(:,1) - q_meas(:,1)).^2));
        sim_metrics.rmse_q2 = sqrt(mean((q_sim(:,2) - q_meas(:,2)).^2));
        sim_metrics.mae_q1 = mean(abs(q_sim(:,1) - q_meas(:,1)));
        sim_metrics.mae_q2 = mean(abs(q_sim(:,2) - q_meas(:,2)));
        if nargin >= 5 && ~isempty(dq_meas)
            sim_metrics.rmse_dq1 = sqrt(mean((dq_sim(:,1) - dq_meas(:,1)).^2));
            sim_metrics.rmse_dq2 = sqrt(mean((dq_sim(:,2) - dq_meas(:,2)).^2));
        end
    end
end

function dxdt = furuta_deriv(x, u, th, r)
    q1 = x(1); q2 = x(2); dq1 = x(3); dq2 = x(4);
    
    % Inertia matrix D(q)
    D = [th(1) + th(2) * (sin(q2)^2),   th(3) * cos(q2);
         th(3) * cos(q2),              th(4)];
     
    % Coriolis and centrifugal vector C(q, dq)*dq
    C_dq = [th(2) * sin(2*q2) * dq1 * dq2 - th(3) * sin(q2) * (dq2^2);
            -0.5 * th(2) * sin(2*q2) * (dq1^2)];
        
    % Gravity vector G(q)
    G = [0; -th(5) * sin(q2)];
    
    % Friction vector F(dq)
    F = [th(6) * dq1 + th(8) * tanh(r * dq1);
         th(7) * dq2 + th(9) * tanh(r * dq2)];
     
    tau = [u; 0];
    
    % System acceleration ddq = D \ (tau - C_dq - G - F)
    ddq = D \ (tau - C_dq - G - F);
    
    dxdt = [dq1; dq2; ddq(1); ddq(2)];
end
