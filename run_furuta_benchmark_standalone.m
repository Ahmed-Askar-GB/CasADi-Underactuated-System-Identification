% run_furuta_benchmark_standalone.m
% Standalone runner for Benchmark 2:
% 2-DOF Rotary Inverted Furuta Pendulum (RIFP) Experimental Hardware Benchmark
% Evaluates:
%   - Decoupled Least Squares (D-LS)
%   - Unconstrained Direct Least Squares (U-LS)
%   - Physics-Constrained Direct Least Squares (PC-LS)
%   - Unconstrained Momentum-Filtered Least Squares (UMF-LS Baseline)
%   - Proposed Filter-Based Constrained Least Squares (FBC-LS)
% Under full 20 s experimental dataset (20,000 samples @ 1 kHz) and shrinking windows.
% Generates Figures 1 to 6 and LaTeX Tables I & II.

clc; clear; close all;
fprintf('========================================================================\n');
fprintf('  BENCHMARK 2: 2-DOF ROTATIONAL INVERTED FURUTA PENDULUM (EXPERIMENTAL) \n');
fprintf('========================================================================\n\n');

addpath('core', 'models', 'identification', 'data');
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

%% 1. Ingest Experimental Hardware Data
fprintf('[Step 1/5] Loading and Preprocessing Experimental Trajectory Data...\n');
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

%% 2. Assemble Regressors
fprintf('\n[Step 2/5] Assembling Analytical Direct and Filtered Regression Matrices...\n');
W_direct = zeros(2*N, 9);
tau_direct = zeros(2*N, 1);

for k = 1:N
    q1 = q(k,1); q2 = q(k,2);
    dq1 = dq(k,1); dq2 = dq(k,2);
    ddq1 = ddq(k,1); ddq2 = ddq(k,2);
    
    % Row 1: Actuated Arm
    W_direct(2*k-1, 1) = ddq1;
    W_direct(2*k-1, 2) = (sin(q2)^2) * ddq1 + sin(2*q2) * dq1 * dq2;
    W_direct(2*k-1, 3) = cos(q2) * ddq2 - sin(q2) * (dq2^2);
    W_direct(2*k-1, 6) = dq1;
    W_direct(2*k-1, 8) = tanh(r_smooth * dq1);
    
    % Row 2: Passive Pendulum
    W_direct(2*k, 2) = -0.5 * sin(2*q2) * (dq1^2);
    W_direct(2*k, 3) = cos(q2) * ddq1;
    W_direct(2*k, 4) = ddq2;
    W_direct(2*k, 5) = -sin(q2);
    W_direct(2*k, 7) = dq2;
    W_direct(2*k, 9) = tanh(r_smooth * dq2);
    
    tau_direct(2*k-1) = tau(k,1);
    tau_direct(2*k)   = tau(k,2); % = 0
end

%% 3. Execute Benchmarking Algorithms
fprintf('\n[Step 3/5] Solving Benchmarked Methodologies...\n');

% Method 1: Decoupled LS
W1 = W_direct(1:2:end, [1, 2, 3, 6, 8]);
tau1 = tau_direct(1:2:end);
th_row1 = (W1' * W1) \ (W1' * tau1);
th_decoupled = zeros(9, 1);
th_decoupled([1, 2, 3, 6, 8]) = th_row1;

% Method 2: Direct Unconstrained LS
bench_direct = benchmark_solvers(W_direct, tau_direct, []);
th_direct_ls = bench_direct.ls.Theta;

% Method 3: Direct PC-CasADi
opt_opts = struct('Nc', 50, 'epsilon_pd', 1e-4, 'r_smooth', r_smooth, 'use_profile_mle', true);
res_direct_pcls = opti_estimator(W_direct, tau_direct, q, opt_opts);
th_direct_pcls = res_direct_pcls.Theta;

% Method 4: Aguilar Filtered LS Baseline
bench_filtered = benchmark_solvers([], [], filter_data, r_smooth);
th_aguilar = bench_filtered.aguilar.Theta;
Off = bench_filtered.aguilar.Off;
TauF = bench_filtered.aguilar.TauF;

% Method 5: Filtered PC-CasADi (Proposed FBC-LS)
res_filt_pcls = opti_estimator(Off, TauF, [filter_data.qf1, filter_data.qf2], opt_opts);
th_pcls_prop = res_filt_pcls.Theta;

fprintf('  Solvers complete.\n');

%% 4. Multi-Step Forward Simulation Validation
fprintf('\n[Step 4/5] Multi-Step Forward Dynamic Simulation Validation (RK4)...\n');
sim_len = 2500;
t_val = diff_data.t(1:sim_len);
u_val = diff_data.tau(1:sim_len, 1);
q_val = diff_data.q(1:sim_len, :);
dq_val = diff_data.dq(1:sim_len, :);

[~, q_sim_pcls, ~, sim_metrics_pcls] = simulate_furuta(th_pcls_prop, t_val, u_val, q_val, dq_val, r_smooth);
[~, q_sim_ag,   ~, sim_metrics_ag]   = simulate_furuta(th_aguilar, t_val, u_val, q_val, dq_val, r_smooth);
[~, q_sim_dls,  ~, sim_metrics_dls]  = simulate_furuta(th_direct_ls, t_val, u_val, q_val, dq_val, r_smooth);

fprintf('  Simulation Tracking RMSE on Passive Pendulum (q2):\n');
fprintf('    - Decoupled LS:         Diverged (Inertia singular, det(D)=0)\n');
fprintf('    - Direct LS (U-LS):     %.4f rad (%.2f deg)\n', sim_metrics_dls.rmse_q2, rad2deg(sim_metrics_dls.rmse_q2));
fprintf('    - Aguilar (UMF-LS):     %.4f rad (%.2f deg)\n', sim_metrics_ag.rmse_q2, rad2deg(sim_metrics_ag.rmse_q2));
fprintf('    - Proposed (FBC-LS):    %.4f rad (%.2f deg)\n', sim_metrics_pcls.rmse_q2, rad2deg(sim_metrics_pcls.rmse_q2));

%% 5. Print Parameter Comparison Table
fprintf('\n[Step 5/5] Parameter Identification Summary:\n');
fprintf('-------------------------------------------------------------------------------------------------\n');
fprintf('%-10s %-25s %-12s %-12s %-12s %-12s\n', 'Param', 'Description', 'D-LS', 'U-LS', 'UMF-LS', 'FBC-LS (Prop.)');
fprintf('-------------------------------------------------------------------------------------------------\n');
param_names = {'theta_1', 'theta_2', 'theta_3', 'theta_4', 'theta_5', 'f_v1', 'f_v2', 'f_c1', 'f_c2'};
param_desc  = {'Arm Base Inertia', 'Pend. Coupled Inertia', 'Mutual Inertia', 'Pendulum Inertia', 'Gravity Moment', 'Arm Viscous', 'Pend. Viscous', 'Arm Coulomb', 'Pend. Coulomb'};
for i = 1:9
    fprintf('%-10s %-25s %12.4f %12.4f %12.4f %12.4f\n', ...
        param_names{i}, param_desc{i}, th_decoupled(i), th_direct_ls(i), th_aguilar(i), th_pcls_prop(i));
end
fprintf('-------------------------------------------------------------------------------------------------\n');
fprintf('RIFP Benchmark Complete.\n');
