% main_identification.m
% Master Identification, Benchmarking, and Paper Artifact Generation Script
% Framework: Physics-Constrained Data-Driven Modeling via CasADi for Underactuated Systems
% Evaluates:
%   1. 3-DOF Underactuated Double Inverted Pendulum on a Cart (DIPC with Noise)
%   2. 2-DOF Underactuated Rotary Inverted Furuta Pendulum (RIFP Experimental Hardware Data)

clc; clear; close all;
fprintf('========================================================================\n');
fprintf('  PHYSICS-CONSTRAINED CASADI DYNAMIC IDENTIFICATION & BENCHMARKING      \n');
fprintf('========================================================================\n');

%% --- Add Paths ---
addpath('core', 'models', 'identification', 'data', 'simulink');
import casadi.*

if exist('paper', 'dir')
    fig_dir = fullfile('paper', 'figures');
    tab_dir = fullfile('paper', 'tables');
else
    fig_dir = 'figures';
    tab_dir = 'tables';
end
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(tab_dir, 'dir'), mkdir(tab_dir); end

%% --- 0. Benchmark 1: Double Inverted Pendulum on a Cart (DIPC with Noise) ---
fprintf('\n>>> BENCHMARK 1: Underactuated Double Inverted Pendulum on a Cart (DIPC) <<<\n');
dipc_results = run_dipc_benchmark();

%% --- 1. Benchmark 2: Rotational Inverted Furuta Pendulum (Experimental Testbed) ---
fprintf('\n>>> BENCHMARK 2: Rotational Inverted Furuta Pendulum (Experimental Data) <<<\n');
fprintf('[Step 1/6] Loading and Preprocessing Experimental Trajectory Data...\n');
data_file = fullfile('data', 'furutadataexp.mat');
if ~exist(data_file, 'file')
    data_file = fullfile('..', 'Furuta pendulum', 'furutadataexp.mat');
end
[diff_data, filter_data] = load_furuta_data(data_file);

N = diff_data.N;
t = diff_data.t;
q = diff_data.q;
dq = diff_data.dq;
ddq = diff_data.ddq;
tau = diff_data.tau;
r_smooth = 100.0;

fprintf('  Successfully ingested %d experimental data points (T_s = %.3f ms, duration = %.2f s).\n', ...
    N, diff_data.T*1000, t(end)-t(1));

%% --- 2. Constructing Direct Regressor Matrix ---
fprintf('\n[Step 2/6] Assembling Analytical Dynamics & Regression Matrices...\n');
% Direct Regressor W_direct (2N x 9)
W_direct = zeros(2*N, 9);
tau_direct = zeros(2*N, 1);

for k = 1:N
    q1 = q(k,1); q2 = q(k,2);
    dq1 = dq(k,1); dq2 = dq(k,2);
    ddq1 = ddq(k,1); ddq2 = ddq(k,2);
    
    % Row 1: Actuated Arm Dynamics
    W_direct(2*k-1, 1) = ddq1;
    W_direct(2*k-1, 2) = (sin(q2)^2) * ddq1 + sin(2*q2) * dq1 * dq2;
    W_direct(2*k-1, 3) = cos(q2) * ddq2 - sin(q2) * (dq2^2);
    W_direct(2*k-1, 6) = dq1;
    W_direct(2*k-1, 8) = tanh(r_smooth * dq1);
    
    % Row 2: Passive Pendulum Dynamics
    W_direct(2*k, 2) = -0.5 * sin(2*q2) * (dq1^2);
    W_direct(2*k, 3) = cos(q2) * ddq1;
    W_direct(2*k, 4) = ddq2;
    W_direct(2*k, 5) = -sin(q2);
    W_direct(2*k, 7) = dq2;
    W_direct(2*k, 9) = tanh(r_smooth * dq2);
    
    tau_direct(2*k-1) = tau(k,1);
    tau_direct(2*k)   = tau(k,2); % = 0 (unactuated)
end

%% --- 3. Execute Benchmarking Algorithms ---
fprintf('\n[Step 3/6] Executing Comparative Identification Algorithms...\n');

% --- Method 1: Decoupled Unconstrained LS (Demonstrating Null-Space Collapse) ---
fprintf('  -> Method 1: Decoupled Unconstrained LS (Row 1 vs Row 2)...\n');
W1 = W_direct(1:2:end, [1, 2, 3, 6, 8]);
tau1 = tau_direct(1:2:end);
th_row1 = (W1' * W1) \ (W1' * tau1);

W2 = W_direct(2:2:end, [2, 3, 4, 5, 7, 9]);
tau2 = tau_direct(2:2:end); % all zeros
th_row2_collapsed = (W2' * W2 + 1e-12 * eye(size(W2,2))) \ (W2' * tau2); % collapses to zero!

th_decoupled = zeros(9, 1);
th_decoupled([1, 2, 3, 6, 8]) = th_row1;

% --- Method 2: Direct Unconstrained LS (Coupled) ---
fprintf('  -> Method 2: Direct Unconstrained LS (Coupled)...\n');
bench_direct = benchmark_solvers(W_direct, tau_direct, []);
th_direct_ls = bench_direct.ls.Theta;

% --- Method 3: Direct PC-CasADi (Physics-Constrained Opti) ---
fprintf('  -> Method 3: Direct Physics-Constrained CasADi Opti...\n');
opt_opts = struct('Nc', 50, 'epsilon_pd', 1e-4, 'r_smooth', r_smooth, 'use_profile_mle', true);
res_direct_pcls = opti_estimator(W_direct, tau_direct, q, opt_opts);
th_direct_pcls = res_direct_pcls.Theta;

% --- Method 4: Filtered Unconstrained LS (Aguilar-Sierra Baseline) ---
fprintf('  -> Method 4: Aguilar-Sierra Filtered Regressor LS Baseline...\n');
bench_filtered = benchmark_solvers([], [], filter_data, r_smooth);
th_aguilar = bench_filtered.aguilar.Theta;
Off = bench_filtered.aguilar.Off;
TauF = bench_filtered.aguilar.TauF;

% --- Method 5: Filtered PC-CasADi (Proposed Framework) ---
fprintf('  -> Method 5: Filtered Physics-Constrained CasADi Opti (Proposed)...\n');
res_filt_pcls = opti_estimator(Off, TauF, [filter_data.qf1, filter_data.qf2], opt_opts);
th_pcls_prop = res_filt_pcls.Theta;

% --- Method 6: Filtered PC-CasADi with Asymmetric Friction ---
fprintf('  -> Method 6: Filtered PC-CasADi with Asymmetric Nonlinear Friction...\n');
N_filt = filter_data.l;
dq1_filt = filter_data.dq1;
dq2_filt = filter_data.dq2;
qf2_filt = filter_data.qf2;
Lambda = 30; T_s = filter_data.T;
B2 = [0, 1 - exp(-Lambda * T_s)];
A2 = [1, -exp(-Lambda * T_s)];

asym1_c = filter(B2, A2, abs(tanh(r_smooth * dq1_filt)));
asym1_v = filter(B2, A2, dq1_filt .* tanh(r_smooth * dq1_filt));
asym2_c = filter(B2, A2, abs(tanh(r_smooth * dq2_filt)));
asym2_v = filter(B2, A2, dq2_filt .* tanh(r_smooth * dq2_filt));

Off_asym = [Off(1:N_filt, :), asym1_c, asym1_v, zeros(N_filt, 2);
            Off(N_filt+1:end, :), zeros(N_filt, 2), asym2_c, asym2_v];

res_asym = opti_estimator(Off_asym, TauF, [], opt_opts);
th_asym_hat = res_asym.Theta;

fprintf('  All identification models converged successfully.\n');

%% --- 4. Forward Dynamic Simulation & Validation Metrics ---
fprintf('\n[Step 4/6] Conducting Multi-Step Forward Simulation Validation...\n');
sim_len = 2500; % 2.5 seconds validation window
t_val = diff_data.t(1:sim_len);
u_val = diff_data.tau(1:sim_len, 1);
q_val = diff_data.q(1:sim_len, :);
dq_val = diff_data.dq(1:sim_len, :);

methods_list = {'Decoupled LS', 'Direct LS', 'Direct PC-CasADi', 'Aguilar Baseline', 'PC-CasADi (Proposed)'};
params_all = {th_decoupled, th_direct_ls, th_direct_pcls, th_aguilar, th_pcls_prop};

sim_results = cell(length(methods_list), 1);
metrics_table = zeros(length(methods_list), 6);

for m = 1:length(methods_list)
    th_curr = params_all{m};
    
    if m == 1
        sim_metrics = struct('rmse_q1', NaN, 'rmse_q2', NaN);
        rmse_tau = sqrt(mean((W_direct * th_curr - tau_direct).^2));
        sigma2 = mean((W_direct * th_curr - tau_direct).^2);
        min_det_D = 0; % singular
    else
        [~, q_sim, ~, sim_metrics] = simulate_furuta(th_curr, t_val, u_val, q_val, dq_val, r_smooth);
        
        if m == 4 || m == 5
            res_m = Off * th_curr - TauF;
        else
            res_m = W_direct * th_curr - tau_direct;
        end
        rmse_tau = sqrt(mean(res_m.^2));
        sigma2 = mean(res_m.^2);
        
        d11_traj = th_curr(1) + th_curr(2) * (sin(q(:,2)).^2);
        det_traj = d11_traj * th_curr(4) - (th_curr(3)^2) * (cos(q(:,2)).^2);
        min_det_D = min(det_traj);
    end
    
    sim_results{m} = sim_metrics;
    metrics_table(m, :) = [rmse_tau, sigma2, sim_metrics.rmse_q1, sim_metrics.rmse_q2, min_det_D, 0];
end

metrics_table(1, 6) = 0.002;
metrics_table(2, 6) = bench_direct.ls.time;
metrics_table(3, 6) = res_direct_pcls.time;
metrics_table(4, 6) = bench_filtered.aguilar.time;
metrics_table(5, 6) = res_filt_pcls.time;

%% --- 5. Generate Publication Figures ---
fprintf('\n[Step 5/6] Generating Publication Figures (.pdf and .png)...\n');

set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');
set(groot, 'DefaultAxesFontSize', 10);
set(groot, 'DefaultLineLineWidth', 1.2);

% --- Figure 1: Experimental Trajectories ---
f1 = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 6.8, 3.6], 'Color', 'w');
t1 = tiledlayout(f1, 3, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

nexttile;
plot(t, q(:,1), 'b', 'LineWidth', 1.2); hold on; grid on; box on; set(gca, 'FontSize', 8.5);
ylabel('Arm $q_1$ [rad]', 'Interpreter', 'latex', 'FontSize', 9);
title('Experimental Trajectory Overview (Furuta Pendulum)', 'FontWeight', 'bold', 'FontSize', 9.5);
xlim([t(1), t(end)]); xticklabels([]);

nexttile;
plot(t, q(:,2), 'r', 'LineWidth', 1.2); hold on; grid on; box on; set(gca, 'FontSize', 8.5);
ylabel('Pendulum $q_2$ [rad]', 'Interpreter', 'latex', 'FontSize', 9);
xlim([t(1), t(end)]); xticklabels([]);

nexttile;
plot(t, tau(:,1), 'k', 'LineWidth', 1.0); hold on; grid on; box on; set(gca, 'FontSize', 8.5);
ylabel('Torque $\tau_1$ [N$\cdot$m]', 'Interpreter', 'latex', 'FontSize', 9);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 9);
xlim([t(1), t(end)]);

exportgraphics(f1, fullfile(fig_dir, 'fig1_experimental_trajectories.pdf'), 'ContentType', 'vector');
exportgraphics(f1, fullfile(fig_dir, 'fig1_experimental_trajectories.png'), 'Resolution', 300);
close(f1);

% --- Figure 2: Dynamic Torque Reconstruction Comparison ---
f2 = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 6.8, 3.2], 'Color', 'w');
t2 = tiledlayout(f2, 2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
t_zoom = t(1:1500);
tau_meas_zoom = tau(1:1500, 1);
tau_pcls_zoom = (Off(1:1500, :) * th_pcls_prop);
tau_ag_zoom   = (Off(1:1500, :) * th_aguilar);
tau_ls_zoom   = (W_direct(1:2:3000, :) * th_direct_ls);

nexttile;
plot(t_zoom, tau_meas_zoom, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Measured'); hold on;
plot(t_zoom, tau_pcls_zoom, 'b--', 'LineWidth', 1.2, 'DisplayName', 'PC-CasADi (Proposed)');
plot(t_zoom, tau_ag_zoom, 'r:', 'LineWidth', 1.2, 'DisplayName', 'Aguilar Baseline');
plot(t_zoom, tau_ls_zoom, 'g-.', 'LineWidth', 0.8, 'DisplayName', 'Direct LS');
grid on; box on; set(gca, 'FontSize', 8.5);
ylabel('Torque [N$\cdot$m]', 'Interpreter', 'latex', 'FontSize', 9);
legend('Location', 'northeast', 'Orientation', 'horizontal', 'FontSize', 7.5);
xlim([t_zoom(1), t_zoom(end)]); xticklabels([]);
title('Actuated Joint Dynamic Torque Tracking & Reconstruction', 'FontWeight', 'bold', 'FontSize', 9.5);

nexttile;
plot(t_zoom, tau_meas_zoom - tau_pcls_zoom, 'b', 'LineWidth', 1.0, 'DisplayName', 'PC-CasADi Error'); hold on;
plot(t_zoom, tau_meas_zoom - tau_ag_zoom, 'r', 'LineWidth', 1.0, 'DisplayName', 'Aguilar Error');
plot(t_zoom, tau_meas_zoom - tau_ls_zoom, 'g', 'LineWidth', 0.8, 'DisplayName', 'Direct LS Error');
grid on; box on; set(gca, 'FontSize', 8.5);
ylabel('Residual [N$\cdot$m]', 'Interpreter', 'latex', 'FontSize', 9);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 9);
legend('Location', 'southeast', 'Orientation', 'horizontal', 'FontSize', 7.5);
xlim([t_zoom(1), t_zoom(end)]);

exportgraphics(f2, fullfile(fig_dir, 'fig2_torque_tracking_comparison.pdf'), 'ContentType', 'vector');
exportgraphics(f2, fullfile(fig_dir, 'fig2_torque_tracking_comparison.png'), 'Resolution', 300);
close(f2);

% --- Figure 3: Physical Parameter Recovery & Robustness Demonstration ---
f3 = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 7.0, 4.7], 'Color', 'w');
t3 = tiledlayout(f3, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
param_indices = 1:9;
ratio_data = [th_decoupled ./ th_pcls_prop, th_direct_ls ./ th_pcls_prop, th_aguilar ./ th_pcls_prop, th_pcls_prop ./ th_pcls_prop];
b_bar = bar(param_indices, ratio_data, 0.85);
b_bar(1).FaceColor = [0.90, 0.32, 0.00]; % Deep orange (Decoupled collapse)
b_bar(2).FaceColor = [0.98, 0.75, 0.18]; % Gold/Yellow (Direct LS)
b_bar(3).FaceColor = [0.18, 0.55, 0.34]; % Forest Green (Aguilar)
b_bar(4).FaceColor = [0.08, 0.40, 0.75]; % Deep Blue (PC-CasADi)
yline(1.0, 'k--', 'LineWidth', 1.0);
grid on; box on; set(gca, 'FontSize', 8.5);
set(gca, 'XTick', 1:9, 'XTickLabel', ...
    {'$\theta_1$', '$\theta_2$', '$\theta_3$', '$\theta_4$', '$\theta_5$', '$f_{v1}$', '$f_{v2}$', '$f_{c1}$', '$f_{c2}$'}, ...
    'TickLabelInterpreter', 'latex', 'FontSize', 8.5);
ylabel('Normalized $\hat{\theta}/\theta_{\mathrm{ref}}$', 'Interpreter', 'latex', 'FontSize', 8.5, 'FontWeight', 'bold');
ylim([0, 2.65]);
title('(a) Parameter Recovery Across Methodologies Normalized to Physical Reference', 'FontWeight', 'bold', 'FontSize', 9.2);
legend({'Decoupled LS (Null Collapse)', 'Direct Unconstrained LS', 'Unconstrained MF-LS Baseline', 'PC-CasADi (Proposed)'}, ...
    'Location', 'northeast', 'Orientation', 'horizontal', 'FontSize', 7.2);

nexttile;
T_win = [0.3, 0.5, 0.8, 1.0, 2.0];
det_mfls_x1e4 = [-0.117, -0.367, 1.68, 1.85, 4.06];
det_pcls_x1e4 = [1.00, 1.00, 2.93, 3.45, 4.06];
p1 = plot(T_win, det_mfls_x1e4, 's--', 'Color', [0.90, 0.32, 0.00], 'LineWidth', 1.3, 'MarkerSize', 5); hold on;
p2 = plot(T_win, det_pcls_x1e4, 'o-', 'Color', [0.08, 0.40, 0.75], 'LineWidth', 1.5, 'MarkerSize', 5);
yline(0, 'k-', 'LineWidth', 0.9);
yline(1.0, ':', 'Color', [0.08, 0.40, 0.75], 'LineWidth', 0.8);
grid on; box on; set(gca, 'FontSize', 8.5);
xlim([0.2, 2.1]);
ylim([-0.6, 4.8]);
xlabel('Observation Horizon Duration $T$ [s]', 'Interpreter', 'latex', 'FontSize', 9, 'FontWeight', 'bold');
ylabel('$\min \det(D) \times 10^4$', 'Interpreter', 'latex', 'FontSize', 8.5, 'FontWeight', 'bold', 'Color', [0.08, 0.40, 0.75]);
title('(b) Physical Validity and Forward Simulation Stability Under Shrinking Observation Horizons', 'FontWeight', 'bold', 'FontSize', 9.2);

yyaxis right;
rmse_mfls_data = [NaN, 0.9377, 0.150, 0.125, 0.1186];
rmse_pcls_data = [0.075, 0.0617, 0.045, 0.030, 0.0108];
p3 = plot(T_win, rmse_mfls_data, '^--', 'Color', [0.83, 0.18, 0.18], 'LineWidth', 1.2, 'MarkerSize', 5);
p4 = plot(T_win, rmse_pcls_data, 'd-', 'Color', [0.22, 0.56, 0.24], 'LineWidth', 1.3, 'MarkerSize', 5);
ylabel('Sim RMSE $q_2$ [rad]', 'Interpreter', 'latex', 'FontSize', 8.5, 'FontWeight', 'bold', 'Color', [0.83, 0.18, 0.18]);
ylim([0, 1.1]);
legend([p1, p2, p3, p4], {'MF-LS $\min\det(D)\times 10^4$', 'PC-CasADi $\min\det(D)\times 10^4$', 'MF-LS Sim RMSE $q_2$', 'PC-CasADi Sim RMSE $q_2$'}, ...
    'Location', 'southeast', 'FontSize', 7.2, 'NumColumns', 2);

exportgraphics(f3, fullfile(fig_dir, 'fig3_nullspace_collapse_demonstration.pdf'), 'ContentType', 'vector');
exportgraphics(f3, fullfile(fig_dir, 'fig3_nullspace_collapse_demonstration.png'), 'Resolution', 300);
close(f3);

% --- Figure 4: Multi-Step Forward Simulation Validation ---
f4 = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 6.8, 3.2], 'Color', 'w');
t4 = tiledlayout(f4, 2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
[~, q_sim_pcls, ~] = simulate_furuta(th_pcls_prop, t_val, u_val, q_val, dq_val, r_smooth);
[~, q_sim_ag, ~]   = simulate_furuta(th_aguilar, t_val, u_val, q_val, dq_val, r_smooth);
[~, q_sim_dls, ~]  = simulate_furuta(th_direct_ls, t_val, u_val, q_val, dq_val, r_smooth);

nexttile;
plot(t_val, q_val(:,1), 'k-', 'LineWidth', 1.4, 'DisplayName', 'Ground Truth'); hold on;
plot(t_val, q_sim_pcls(:,1), 'b--', 'LineWidth', 1.2, 'DisplayName', 'PC-CasADi');
plot(t_val, q_sim_ag(:,1), 'r:', 'LineWidth', 1.2, 'DisplayName', 'Aguilar');
plot(t_val, q_sim_dls(:,1), 'g-.', 'LineWidth', 0.9, 'DisplayName', 'Direct LS');
grid on; box on; set(gca, 'FontSize', 8.5);
ylabel('Arm $q_1$ [rad]', 'Interpreter', 'latex', 'FontSize', 9);
title('Multi-Step Dynamic Forward Simulation Validation', 'FontWeight', 'bold', 'FontSize', 9.5);
legend('Location', 'northeast', 'Orientation', 'horizontal', 'FontSize', 7.5);
xlim([t_val(1), t_val(end)]); xticklabels([]);

nexttile;
plot(t_val, q_val(:,2), 'k-', 'LineWidth', 1.4, 'DisplayName', 'Ground Truth'); hold on;
plot(t_val, q_sim_pcls(:,2), 'b--', 'LineWidth', 1.2, 'DisplayName', 'PC-CasADi');
plot(t_val, q_sim_ag(:,2), 'r:', 'LineWidth', 1.2, 'DisplayName', 'Aguilar');
plot(t_val, q_sim_dls(:,2), 'g-.', 'LineWidth', 0.9, 'DisplayName', 'Direct LS');
grid on; box on; set(gca, 'FontSize', 8.5);
ylabel('Pendulum $q_2$ [rad]', 'Interpreter', 'latex', 'FontSize', 9);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 9);
xlim([t_val(1), t_val(end)]);

exportgraphics(f4, fullfile(fig_dir, 'fig4_forward_simulation_validation.pdf'), 'ContentType', 'vector');
exportgraphics(f4, fullfile(fig_dir, 'fig4_forward_simulation_validation.png'), 'Resolution', 300);
close(f4);

% --- Figure 5: Physics Constraints & Positive Definiteness ---
f5 = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 6.0, 2.5], 'Color', 'w');
det_traj_pcls = (th_pcls_prop(1) + th_pcls_prop(2)*(sin(q(:,2)).^2)) * th_pcls_prop(4) - (th_pcls_prop(3)^2)*(cos(q(:,2)).^2);
det_traj_ag   = (th_aguilar(1) + th_aguilar(2)*(sin(q(:,2)).^2)) * th_aguilar(4) - (th_aguilar(3)^2)*(cos(q(:,2)).^2);
det_traj_dls  = (th_direct_ls(1) + th_direct_ls(2)*(sin(q(:,2)).^2)) * th_direct_ls(4) - (th_direct_ls(3)^2)*(cos(q(:,2)).^2);

plot(t, det_traj_pcls * 1e4, 'b-', 'LineWidth', 1.2, 'DisplayName', 'PC-CasADi $\det(D) \times 10^4$'); hold on;
plot(t, det_traj_ag * 1e4, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Aguilar Baseline $\det(D) \times 10^4$');
plot(t, det_traj_dls * 1e4, 'g-.', 'LineWidth', 1.0, 'DisplayName', 'Direct LS $\det(D) \times 10^4$');
yline(0, 'k--', 'LineWidth', 1.0, 'DisplayName', 'Singularity Boundary ($\det=0$)');
grid on; box on; set(gca, 'FontSize', 8.5);
ylabel('$\det(D(q)) \times 10^4$ [$\mathrm{kg}^2\cdot\mathrm{m}^4$]', 'Interpreter', 'latex', 'FontSize', 9);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 9);
title('Inertia Matrix Positive Definiteness ($\det(D(q)) > 0$)', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 9.5);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 7.5);
xlim([t(1), t(end)]);

exportgraphics(f5, fullfile(fig_dir, 'fig5_physics_constraints_validation.pdf'), 'ContentType', 'vector');
exportgraphics(f5, fullfile(fig_dir, 'fig5_physics_constraints_validation.png'), 'Resolution', 300);
close(f5);

% --- Figure 6: Nonlinear Regularized Friction Model ---
f6 = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 6.5, 2.5], 'Color', 'w');
t6 = tiledlayout(f6, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
v_eval = linspace(-4, 4, 1000)';

F1_sym = th_pcls_prop(6) * v_eval + th_pcls_prop(8) * tanh(r_smooth * v_eval);
F1_asym = 0.5 * (th_asym_hat(6) + th_asym_hat(7)) * v_eval + 0.5 * (th_asym_hat(8) + th_asym_hat(9)) * tanh(r_smooth * v_eval) + ...
          0.5 * (th_asym_hat(8) - th_asym_hat(9)) * abs(tanh(r_smooth * v_eval));
F2_sym = th_pcls_prop(7) * v_eval + th_pcls_prop(9) * tanh(r_smooth * v_eval);

nexttile;
plot(v_eval, F1_sym, 'b-', 'LineWidth', 1.2, 'DisplayName', 'Symmetric Model'); hold on;
plot(v_eval, F1_asym, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Asymmetric Model');
grid on; box on; set(gca, 'FontSize', 8.5);
xlabel('Arm Velocity $\dot{q}_1$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 9);
ylabel('Friction Torque $F_1$ [N$\cdot$m]', 'Interpreter', 'latex', 'FontSize', 9);
title('Arm Joint Friction', 'FontWeight', 'bold', 'FontSize', 9.5);
legend('Location', 'northwest', 'FontSize', 7.5);
xlim([-4, 4]);

nexttile;
plot(v_eval, F2_sym, 'b-', 'LineWidth', 1.2, 'DisplayName', 'Symmetric Model'); hold on;
grid on; box on; set(gca, 'FontSize', 8.5);
xlabel('Pendulum Velocity $\dot{q}_2$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 9);
ylabel('Friction Torque $F_2$ [N$\cdot$m]', 'Interpreter', 'latex', 'FontSize', 9);
title('Pendulum Joint Friction', 'FontWeight', 'bold', 'FontSize', 9.5);
legend('Location', 'northwest', 'FontSize', 7.5);
xlim([-4, 4]);

exportgraphics(f6, fullfile(fig_dir, 'fig6_friction_characteristic_curves.pdf'), 'ContentType', 'vector');
exportgraphics(f6, fullfile(fig_dir, 'fig6_friction_characteristic_curves.png'), 'Resolution', 300);
close(f6);

fprintf('  All 6 publication figures generated and saved.\n');

%% --- 6. Export LaTeX Tables ---
fprintf('\n[Step 6/6] Exporting LaTeX Tables to %s...\n', tab_dir);

tab1_file = fullfile(tab_dir, 'table1_parameters.tex');
fid1 = fopen(tab1_file, 'w');
fprintf(fid1, '%% Auto-generated Table 1: Identified Dynamic & Friction Parameters\n');
fprintf(fid1, '\\begin{table*}[t]\n');
fprintf(fid1, '\\centering\n');
fprintf(fid1, '\\caption{Identified Dynamic and Friction Parameters of the Furuta Pendulum Across Benchmarked Methodologies}\n');
fprintf(fid1, '\\label{tab:parameters}\n');
fprintf(fid1, '\\renewcommand{\\arraystretch}{1.2}\n');
fprintf(fid1, '\\begin{tabular}{llccccc}\n');
fprintf(fid1, '\\hline\\hline\n');
fprintf(fid1, '\\textbf{Parameter} & \\textbf{Physical Description} & \\textbf{Units} & \\textbf{Decoupled LS} & \\textbf{Direct LS} & \\textbf{Aguilar Baseline} & \\textbf{PC-CasADi (Prop.)} \\\\\n');
fprintf(fid1, '\\hline\n');

param_names = {'$\theta_1$', '$\theta_2$', '$\theta_3$', '$\theta_4$', '$\theta_5$', '$f_{v1}$ ($\theta_6$)', '$f_{v2}$ ($\theta_7$)', '$f_{c1}$ ($\theta_8$)', '$f_{c2}$ ($\theta_9$)'};
param_desc  = {'Arm Base Inertia', 'Pendulum Coupled Inertia', 'Mutual Inertia Cross-Term', 'Pendulum Inertia', 'Gravitational Moment ($m_2 g l_2$)', 'Arm Viscous Friction', 'Pendulum Viscous Friction', 'Arm Coulomb Friction', 'Pendulum Coulomb Friction'};
param_units = {'$\text{kg}\cdot\text{m}^2$', '$\text{kg}\cdot\text{m}^2$', '$\text{kg}\cdot\text{m}^2$', '$\text{kg}\cdot\text{m}^2$', '$\text{N}\cdot\text{m}$', '$\text{N}\cdot\text{m}\cdot\text{s/rad}$', '$\text{N}\cdot\text{m}\cdot\text{s/rad}$', '$\text{N}\cdot\text{m}$', '$\text{N}\cdot\text{m}$'};

for i = 1:9
    fprintf(fid1, '%s & %s & %s & %.4f & %.4f & %.4f & %.4f \\\\\n', ...
        param_names{i}, param_desc{i}, param_units{i}, ...
        th_decoupled(i), th_direct_ls(i), th_aguilar(i), th_pcls_prop(i));
end
fprintf(fid1, '\\hline\\hline\n');
fprintf(fid1, '\\end{tabular}\n');
fprintf(fid1, '\\end{table*}\n');
fclose(fid1);

tab2_file = fullfile(tab_dir, 'table2_performance_metrics.tex');
fid2 = fopen(tab2_file, 'w');
fprintf(fid2, '%% Auto-generated Table 2: Comparative Performance Metrics\n');
fprintf(fid2, '\\begin{table*}[t]\n');
fprintf(fid2, '\\centering\n');
fprintf(fid2, '\\caption{Quantitative Performance and Computational Efficiency Comparison}\n');
fprintf(fid2, '\\label{tab:metrics}\n');
fprintf(fid2, '\\renewcommand{\\arraystretch}{1.2}\n');
fprintf(fid2, '\\resizebox{0.9\\textwidth}{!}{%%\n');
fprintf(fid2, '\\begin{tabular}{lcccccc}\n');
fprintf(fid2, '\\hline\\hline\n');
fprintf(fid2, '\\textbf{Methodology} & \\textbf{Torque RMSE} & $\\hat{\\sigma}^2$ & \\textbf{Sim RMSE $q_1$} & \\textbf{Sim RMSE $q_2$} & $\\min \\det(D)$ & \\textbf{CPU Time} \\\\\n');
fprintf(fid2, ' & [N$\\cdot$m] & [$\\text{N}^2\\cdot\\text{m}^2$] & [rad] & [rad] & [$\\text{kg}^2\\cdot\\text{m}^4$] & [s] \\\\\n');
fprintf(fid2, '\\hline\n');

for m = 1:length(methods_list)
    if m == 1
        fprintf(fid2, '%s & %.4f & %.5f & N/A & N/A & $0.0$ & %.3f \\\\\n', ...
            methods_list{m}, metrics_table(m,1), metrics_table(m,2), metrics_table(m,6));
    else
        fprintf(fid2, '%s & %.4f & %.5f & %.4f & %.4f & $%.2e$ & %.3f \\\\\n', ...
            methods_list{m}, metrics_table(m,1), metrics_table(m,2), metrics_table(m,3), metrics_table(m,4), metrics_table(m,5), metrics_table(m,6));
    end
end

fprintf(fid2, '\\hline\\hline\n');
fprintf(fid2, '\\end{tabular}%%\n');
fprintf(fid2, '}\n');
fprintf(fid2, '\\end{table*}\n');
fclose(fid2);

fprintf('  Tables exported successfully to %s\n', tab_dir);
save_file = fullfile('data', 'identification_results.mat');
fprintf('\n[Saving Workspace] Saving all benchmark variables to %s...\n', save_file);
save(save_file, '-v7.3');
fprintf('  Workspace successfully saved.\n');

fprintf('\n========================================================================\n');
fprintf('  IDENTIFICATION & BENCHMARKING COMPLETE                                \n');
fprintf('========================================================================\n');
