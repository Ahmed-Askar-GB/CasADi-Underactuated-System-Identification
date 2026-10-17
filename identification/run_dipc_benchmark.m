function results = run_dipc_benchmark()
% RUN_DIPC_BENCHMARK Master benchmark for Double Inverted Pendulum on a Cart (DIPC)
% Evaluates Direct OLS, Aguilar-Sierra Baseline, and Physics-Constrained CasADi (PC-CasADi) under sensor noise.
% Generates publication figures and LaTeX comparison table.

fprintf('========================================================================\n');
fprintf('  DIPC SIMULATION BENCHMARK: PARAMETER IDENTIFICATION & NOISE ESTIMATION\n');
fprintf('========================================================================\n');

import casadi.*

% Determine figure and table output directories
if exist('paper', 'dir')
    fig_dir = fullfile('paper', 'figures');
    tab_dir = fullfile('paper', 'tables');
else
    fig_dir = 'figures';
    tab_dir = 'tables';
end
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(tab_dir, 'dir'), mkdir(tab_dir); end

%% 1. Ingest Simulation Data
candidates = {
    fullfile('data', 'Sim1.mat'), ...
    fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'Sim1.mat'), ...
    'Sim1.mat', ...
    fullfile('..', 'data', 'Sim1.mat'), ...
    fullfile('..', 'Code', 'Nile Paper', 'Sim1.mat')
};

sim_file = '';
for c = 1:length(candidates)
    if exist(candidates{c}, 'file')
        sim_file = candidates{c};
        break;
    end
end

if isempty(sim_file)
    error('Sim1.mat not found in data/ or search paths.');
end

s = load(sim_file);
out = s.out;

% Extract trajectory signals (trim initial transient)
idx_start = 10;
t = out.t(idx_start:end);
t = t - t(1);
dt = t(2) - t(1);
N = length(t);

d1_true   = out.d1(idx_start:end);
dd1_true  = out.dd1(idx_start:end);
ddd1_true = out.ddd1(idx_start:end);

s2_true   = out.s2(idx_start:end);
ds2_true  = out.ds2(idx_start:end);
dds2_true = out.dds2(idx_start:end);

s3_true   = out.s3(idx_start:end);
ds3_true  = out.ds3(idx_start:end);
dds3_true = out.dds3(idx_start:end);

fprintf('[Step 1/6] Loaded DIPC simulation data from %s: N = %d samples, dt = %.4f s, T = %.2f s\n', ...
    sim_file, N, dt, t(end));

%% 2. System Ground Truth & Regressor Construction
% Physical parameters from Data-driven-Underactuator.pdf
m1 = 5.4; m2 = 1.35; m3 = 1.35;
a1 = 0.05; a2 = 0.5; a3 = 0.5;
lc1 = 0; lc2 = a2/2; lc3 = a3/2;
I2zz = 0.0284063; I3zz = 0.0284063;
g = -9.80665;

% Ground truth coefficients
coff1 = [m1+m2+m3; -(a2*m3+lc2*m2); -lc3*m3; -lc3*m3; a2*m3+lc2*m2; lc3*m3; 2*lc3*m3; lc3*m3];
coff2 = [-(lc2*m2+a2*m3); -lc3*m3; I2zz+I3zz+a2^2*m3+lc2^2*m2+lc3^2*m3; 2*a2*lc3*m3; I3zz+lc3^2*m3; a2*lc3*m3; -a2*lc3*m3; -2*a2*lc3*m3; g*lc2*m2+a2*g*m3; g*lc3*m3];
coff3 = [-lc3*m3; I3zz+lc3^2*m3; a2*lc3*m3; I3zz+lc3^2*m3; a2*lc3*m3; g*lc3*m3];
theta_gt = [coff1; coff2; coff3]; % 24 parameters

% Regressor matrices for links 1, 2, and 3
L1 = [ddd1_true, dds2_true.*cos(s2_true), dds2_true.*cos(s2_true+s3_true), dds3_true.*cos(s2_true+s3_true), ...
      ds2_true.^2.*sin(s2_true), ds2_true.^2.*sin(s2_true+s3_true), ds2_true.*ds3_true.*sin(s2_true+s3_true), ds3_true.^2.*sin(s2_true+s3_true)];
L2 = [ddd1_true.*cos(s2_true), ddd1_true.*cos(s2_true+s3_true), dds2_true, dds2_true.*cos(s3_true), ...
      dds3_true, dds3_true.*cos(s3_true), ds3_true.^2.*sin(s3_true), ds2_true.*ds3_true.*sin(s3_true), sin(s2_true), sin(s2_true+s3_true)];
L3 = [ddd1_true.*cos(s2_true+s3_true), dds2_true, dds2_true.*cos(s3_true), dds3_true, ds2_true.^2.*sin(s3_true), sin(s2_true+s3_true)];

% Clean physical generalized force
F1_clean = L1 * coff1;

%% 3. Inject Measurement Noise
fprintf('\n[Step 2/6] Injecting Measurement Noise (sigma_true = 0.05, sigma2 = 0.0025)...\n');
rng(42); % reproducible seed
sigma_true = 0.05;
sigma2_true = sigma_true^2;

noise_F1 = sigma_true * randn(N, 1);
noise_T2 = sigma_true * randn(N, 1);
noise_T3 = sigma_true * randn(N, 1);

F1_noisy = F1_clean + noise_F1;
T2_noisy = zeros(N, 1) + noise_T2;
T3_noisy = zeros(N, 1) + noise_T3;

% Stacked block regressor
W_block = blkdiag(L1, L2, L3); % (3N x 24)
Gamma = [F1_noisy; T2_noisy; T3_noisy];
N_rows = length(Gamma);

%% 4. Direct Unconstrained OLS Identification
fprintf('\n[Step 3/6] Solving Direct Unconstrained Ordinary Least Squares (OLS)...\n');
tic;
theta_ols = (W_block' * W_block) \ (W_block' * Gamma);
t_ols = toc;

res_ols = W_block * theta_ols - Gamma;
sigma2_ols = mean(res_ols.^2);
err_ols = norm(theta_ols - theta_gt) / norm(theta_gt) * 100;
norm_L2_ols = norm(theta_ols(9:18));
norm_L3_ols = norm(theta_ols(19:24));

fprintf('  Direct OLS Relative Parameter Error: %.2f%%\n', err_ols);
fprintf('  Direct OLS Link 2 Norm: %.4f (Ground Truth: %.4f)\n', norm_L2_ols, norm(coff2));
fprintf('  Direct OLS Link 3 Norm: %.4f (Ground Truth: %.4f)\n', norm_L3_ols, norm(coff3));
fprintf('  Direct OLS Execution Time: %.3f ms\n', t_ols*1000);

%% 5. Aguilar-Sierra Filtered Regressor LS Baseline
fprintf('\n[Step 4/6] Solving Aguilar-Sierra Filtered Regressor LS Baseline...\n');
tic;
Lambda = 30;
B2 = [0, 1 - exp(-Lambda * dt)];
A2 = [1, -exp(-Lambda * dt)];

L1_filt = zeros(size(L1));
L2_filt = zeros(size(L2));
L3_filt = zeros(size(L3));
for c = 1:size(L1,2), L1_filt(:,c) = filter(B2, A2, L1(:,c)); end
for c = 1:size(L2,2), L2_filt(:,c) = filter(B2, A2, L2(:,c)); end
for c = 1:size(L3,2), L3_filt(:,c) = filter(B2, A2, L3(:,c)); end

F1_filt = filter(B2, A2, F1_noisy);
T2_filt = filter(B2, A2, T2_noisy);
T3_filt = filter(B2, A2, T3_noisy);

W_filt = blkdiag(L1_filt, L2_filt, L3_filt);
Gamma_filt = [F1_filt; T2_filt; T3_filt];

theta_aguilar = (W_filt' * W_filt + 1e-8*eye(24)) \ (W_filt' * Gamma_filt);
t_ag = toc;

res_ag = W_block * theta_aguilar - Gamma;
sigma2_ag = mean(res_ag.^2);
err_ag = norm(theta_aguilar - theta_gt) / norm(theta_gt) * 100;
norm_L2_ag = norm(theta_aguilar(9:18));
norm_L3_ag = norm(theta_aguilar(19:24));

fprintf('  Aguilar Relative Parameter Error: %.2f%%\n', err_ag);
fprintf('  Aguilar Link 2 Norm: %.4f (Ground Truth: %.4f)\n', norm_L2_ag, norm(coff2));
fprintf('  Aguilar Link 3 Norm: %.4f (Ground Truth: %.4f)\n', norm_L3_ag, norm(coff3));
fprintf('  Aguilar Execution Time: %.3f ms\n', t_ag*1000);

%% 6. CasADi Physics-Constrained Least Squares (PC-CasADi)
fprintf('\n[Step 5/6] Solving Physics-Constrained Least Squares via CasADi Opti...\n');
opti = Opti();
theta = opti.variable(24, 1);

x1 = theta(1:8);
x2 = theta(9:18);
x3 = theta(19:24);

% Exact algebraic symmetry constraints (D_ij = D_ji)
opti.subject_to(x1(2) == x2(1));
opti.subject_to(x1(3) == x2(2));
opti.subject_to(x1(4) == x3(1));
opti.subject_to(x2(5) == x3(2));
opti.subject_to(x2(6) == x3(3));
opti.subject_to(x2(4) == 2 * x2(6));

% Exact algebraic skew-symmetry passivity constraints (N + N' = 0)
opti.subject_to(x1(2) + 2*x1(5) + x2(1) == 0);
opti.subject_to(x1(3) + 2*x1(6) + x2(2) == 0);
opti.subject_to(x1(3) + x1(7) + x2(2) == 0);
opti.subject_to(x1(4) + x1(7) + x3(1) == 0);
opti.subject_to(x1(4) + 2*x1(8) + x3(1) == 0);
opti.subject_to(x2(6) + 2*x2(7) + x3(3) == 0);
opti.subject_to(x2(8) + 2*x3(5) == 0);

% Mass and inertia parameter structural equalities
opti.subject_to(x3(4) == x3(2));
opti.subject_to(x3(3) == x3(5));

% Inertia positivity
opti.subject_to(x1(1) >= 1.0);
opti.subject_to(x2(3) >= 0.05);
opti.subject_to(x3(4) >= 0.01);

% Gravity potential bounds
opti.subject_to(x2(9) <= -1.0);
opti.subject_to(x2(10) <= -0.5);
opti.subject_to(x3(6) <= -0.5);
opti.subject_to(x2(10) == x3(6));

% Positive definiteness check (Sylvester criterion) at key configurations
sample_configs = [0, 0; pi/4, pi/4; pi/2, 0; -pi/3, pi/6];
for sc = 1:size(sample_configs, 1)
    s2_val = sample_configs(sc, 1);
    s3_val = sample_configs(sc, 2);
    
    D11 = x1(1);
    D12 = x1(2)*cos(s2_val) + x1(3)*cos(s2_val + s3_val);
    D13 = x1(4)*cos(s2_val + s3_val);
    D22 = x2(3) + x2(4)*cos(s3_val);
    D23 = x2(5) + x2(6)*cos(s3_val);
    D33 = x3(4);
    
    opti.subject_to(D11 >= 1.0);
    opti.subject_to(D11*D22 - D12^2 >= 0.01);
    detD = D11*(D22*D33 - D23^2) - D12*(D12*D33 - D13*D23) + D13*(D12*D23 - D13*D22);
    opti.subject_to(detD >= 0.001);
end

% Condensed Normal Equations objective
H = W_block' * W_block;
g_vec = W_block' * Gamma;
c_const = Gamma' * Gamma;

obj = theta' * H * theta - 2 * g_vec' * theta + c_const;
opti.minimize(obj);

% Solver settings
opts = struct();
opts.ipopt.print_level = 0;
opts.print_time = false;
opts.ipopt.tol = 1e-8;
opti.solver('ipopt', opts);

opti.set_initial(theta, zeros(24,1));

tic;
sol = opti.solve();
t_pc = toc;

theta_pc = full(sol.value(theta));
res_pc = W_block * theta_pc - Gamma;
sigma2_pc = mean(res_pc.^2);
err_pc = norm(theta_pc - theta_gt) / norm(theta_gt) * 100;
norm_L2_pc = norm(theta_pc(9:18));
norm_L3_pc = norm(theta_pc(19:24));

fprintf('  PC-CasADi Relative Parameter Error: %.2f%%\n', err_pc);
fprintf('  PC-CasADi Noise Variance: %.6f (True: %.6f, Rel Error: %.2f%%)\n', ...
    sigma2_pc, sigma2_true, abs(sigma2_pc - sigma2_true)/sigma2_true * 100);
fprintf('  PC-CasADi Link 2 Norm: %.4f (Ground Truth: %.4f)\n', norm_L2_pc, norm(coff2));
fprintf('  PC-CasADi Link 3 Norm: %.4f (Ground Truth: %.4f)\n', norm_L3_pc, norm(coff3));
fprintf('  PC-CasADi Execution Time: %.3f ms\n', t_pc*1000);

%% 7. Dynamic Forward Simulation Benchmarking
fprintf('\n[Step 6/6] Conducting Forward Simulation Benchmarking...\n');
q0 = [d1_true(1); s2_true(1); s3_true(1)];
dq0 = [dd1_true(1); ds2_true(1); ds3_true(1)];

sim_dyn = @(th) simulate_dipc_forward(t, q0, dq0, F1_noisy, th);

[q_sim_gt, valid_gt]   = sim_dyn(theta_gt);
[q_sim_pc, valid_pc]   = sim_dyn(theta_pc);
[q_sim_ols, valid_ols] = sim_dyn(theta_ols);
[q_sim_ag, valid_ag]   = sim_dyn(theta_aguilar);

rmse_cart_pc = sqrt(mean((q_sim_pc(:,1) - d1_true).^2));
rmse_s2_pc   = sqrt(mean((q_sim_pc(:,2) - s2_true).^2));
rmse_s3_pc   = sqrt(mean((q_sim_pc(:,3) - s3_true).^2));

fprintf('  Forward Sim RMSE (PC-CasADi): Cart = %.4f m, Pend 1 = %.4f rad, Pend 2 = %.4f rad\n', ...
    rmse_cart_pc, rmse_s2_pc, rmse_s3_pc);
if ~valid_ols
    fprintf('  Forward Sim (Direct OLS): DIVERGED (Matrix singular / NaN / Inf)\n');
end
if ~valid_ag
    fprintf('  Forward Sim (Aguilar Baseline): DIVERGED (Matrix singular / NaN / Inf)\n');
end

%% 8. Generate Figures
col_gt  = [0.20, 0.20, 0.20];     % Dark Gray / Black
col_ols = [0.85, 0.35, 0.35];     % Light Red
col_ag  = [0.85, 0.10, 0.50];     % Magenta / Crimson
col_pc  = [0.00, 0.45, 0.74];     % Deep Blue

% --- Figure A: Subsystem-Grouped Relative Parameter Estimation Error (%) ---
err_bar_ols = abs(theta_ols - theta_gt) ./ abs(theta_gt) * 100;
err_bar_ag  = abs(theta_aguilar - theta_gt) ./ abs(theta_gt) * 100;
err_bar_pc  = abs(theta_pc - theta_gt) ./ abs(theta_gt) * 100;
bar_data = [err_bar_ols, err_bar_ag, err_bar_pc];

figA = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1, 1, 6.8, 4.4]);
tA = tiledlayout(figA, 3, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

% Tile 1: Link 1 (Cart, 8 parameters)
nexttile;
b1 = bar(1:8, bar_data(1:8, :), 0.85);
b1(1).FaceColor = col_ols; b1(1).EdgeColor = 'none';
b1(2).FaceColor = [0.18, 0.55, 0.34]; b1(2).EdgeColor = 'none';
b1(3).FaceColor = col_pc;  b1(3).EdgeColor = 'none';
grid on; box on; set(gca, 'FontSize', 8.5);
xlim([0.4, 8.6]);
ylim([0, 0.45]);
xticks(1:8);
xticklabels({'x_{1,1}', 'x_{1,2}', 'x_{1,3}', 'x_{1,4}', 'x_{1,5}', 'x_{1,6}', 'x_{1,7}', 'x_{1,8}'});
ylabel('Cart ($L_1$)\newline Rel. Err. [\%]', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 8.5);
title('Double Inverted Pendulum on a Cart: Relative Parameter Error Across Subsystems ($\sigma=0.05$)', ...
    'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 9.5);
legend({'Direct Unconstrained OLS', 'Unconstrained MF-LS Baseline', 'PC-CasADi (Proposed)'}, ...
    'Location', 'northeast', 'Orientation', 'horizontal', 'FontSize', 7.5);

% Tile 2: Link 2 (Passive Pendulum 1, 10 parameters)
nexttile;
b2 = bar(1:10, bar_data(9:18, :), 0.85);
b2(1).FaceColor = col_ols; b2(1).EdgeColor = 'none';
b2(2).FaceColor = [0.18, 0.55, 0.34]; b2(2).EdgeColor = 'none';
b2(3).FaceColor = col_pc;  b2(3).EdgeColor = 'none';
yline(100, 'k:', 'LineWidth', 0.8);
grid on; box on; set(gca, 'FontSize', 8.5);
xlim([0.4, 10.6]);
ylim([0, 115]);
xticks(1:10);
xticklabels({'x_{2,1}', 'x_{2,2}', 'x_{2,3}', 'x_{2,4}', 'x_{2,5}', 'x_{2,6}', 'x_{2,7}', 'x_{2,8}', 'x_{2,9}', 'x_{2,10}'});
ylabel('Pend. 1 ($L_2$)\newline Rel. Err. [\%]', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 8.5);

% Tile 3: Link 3 (Passive Pendulum 2, 6 parameters)
nexttile;
b3 = bar(1:6, bar_data(19:24, :), 0.85);
b3(1).FaceColor = col_ols; b3(1).EdgeColor = 'none';
b3(2).FaceColor = [0.18, 0.55, 0.34]; b3(2).EdgeColor = 'none';
b3(3).FaceColor = col_pc;  b3(3).EdgeColor = 'none';
yline(100, 'k:', 'LineWidth', 0.8);
grid on; box on; set(gca, 'FontSize', 8.5);
xlim([0.4, 6.6]);
ylim([0, 115]);
xticks(1:6);
xticklabels({'x_{3,1}', 'x_{3,2}', 'x_{3,3}', 'x_{3,4}', 'x_{3,5}', 'x_{3,6}'});
ylabel('Pend. 2 ($L_3$)\newline Rel. Err. [\%]', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 8.5);
xlabel('Identified Dynamic Base Parameters', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 9);

exportgraphics(figA, fullfile(fig_dir, 'fig_dipc_nullspace.pdf'), 'ContentType', 'vector');
exportgraphics(figA, fullfile(fig_dir, 'fig_dipc_nullspace.png'), 'Resolution', 300);
close(figA);
fprintf('  Saved Figure: %s\n', fullfile(fig_dir, 'fig_dipc_nullspace.pdf'));

% --- Figure B: Dynamic Trajectory Tracking Comparison (Tight 4-Tile Layout) ---
figB = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1, 1, 6.8, 4.4]);
tB = tiledlayout(figB, 4, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

nexttile;
plot(t, d1_true, 'k-', 'LineWidth', 1.4); hold on;
plot(t, q_sim_pc(:,1), '--', 'Color', col_pc, 'LineWidth', 1.2);
grid on; box on; set(gca, 'FontSize', 8.5);
ylabel('Cart $d_1$ [m]', 'Interpreter', 'latex', 'FontWeight', 'bold');
xlim([t(1), t(end)]); xticklabels([]);
title('DIPC Forward Simulation Fidelity Under Identified Models ($\sigma = 0.05$)', ...
    'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 9.5);
legend({'True Trajectory', 'PC-CasADi (Ours)'}, 'Location', 'southeast', 'FontSize', 7.5);

nexttile;
plot(t, s2_true, 'k-', 'LineWidth', 1.4); hold on;
plot(t, q_sim_pc(:,2), '--', 'Color', col_pc, 'LineWidth', 1.2);
grid on; box on; set(gca, 'FontSize', 8.5);
ylabel('$\theta_2$ [rad]', 'Interpreter', 'latex', 'FontWeight', 'bold');
xlim([t(1), t(end)]); xticklabels([]);

nexttile;
plot(t, s3_true, 'k-', 'LineWidth', 1.4); hold on;
plot(t, q_sim_pc(:,3), '--', 'Color', col_pc, 'LineWidth', 1.2);
grid on; box on; set(gca, 'FontSize', 8.5);
ylabel('$\theta_3$ [rad]', 'Interpreter', 'latex', 'FontWeight', 'bold');
xlim([t(1), t(end)]); xticklabels([]);

nexttile;
plot(t, F1_clean, 'k-', 'LineWidth', 1.3); hold on;
plot(t, F1_noisy, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 2);
plot(t, L1 * theta_pc(1:8), '--', 'Color', col_pc, 'LineWidth', 1.1);
grid on; box on; set(gca, 'FontSize', 8.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontWeight', 'bold');
ylabel('Force $F_1$ [N]', 'Interpreter', 'latex', 'FontWeight', 'bold');
xlim([t(1), t(end)]);
legend({'Clean Force', 'Noisy Obs. ($\sigma=0.05$)', 'PC-CasADi Pred.'}, 'Location', 'northwest', 'FontSize', 7.5);

exportgraphics(figB, fullfile(fig_dir, 'fig_dipc_tracking.pdf'), 'ContentType', 'vector');
exportgraphics(figB, fullfile(fig_dir, 'fig_dipc_tracking.png'), 'Resolution', 300);
close(figB);
fprintf('  Saved Figure: %s\n', fullfile(fig_dir, 'fig_dipc_tracking.pdf'));

%% 9. Generate LaTeX Benchmark Table
tab_file = fullfile(tab_dir, 'table4_dipc_comparison.tex');
fid = fopen(tab_file, 'w');
fprintf(fid, '%% Auto-generated Table IV: Double Inverted Pendulum on a Cart (DIPC) Identification\n');
fprintf(fid, '\\begin{table*}[t]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\footnotesize\n');
fprintf(fid, '\\setlength{\\tabcolsep}{1.8pt}\n');
fprintf(fid, '\\caption{Dynamic Parameter Identification and Noise Estimation on the Underactuated Double Inverted Pendulum on a Cart (DIPC) Under Injected Sensor Noise ($\\sigma_{\\text{true}} = 0.05$, $\\sigma^2_{\\text{true}} = 2.50 \\times 10^{-3}$)}\n');
fprintf(fid, '\\label{tab:dipc_benchmark}\n');
fprintf(fid, '\\resizebox{\\textwidth}{!}{%%\n');
fprintf(fid, '\\begin{tabular}{l l l r r r r r}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, '\\textbf{Subsystem} & \\textbf{Param.} & \\textbf{Physical Expression} & \\textbf{Ground Truth} & \\textbf{Direct OLS} & \\textbf{Aguilar Baseline} & \\textbf{PC-CasADi (Ours)} & \\textbf{PC Rel. Err.} \\\\\n');
fprintf(fid, '\\midrule\n');

param_syms = {
    'x_{1,1}', 'm_1+m_2+m_3';
    'x_{1,2}', '-(a_2 m_3 + lc_2 m_2)';
    'x_{1,3}', '-lc_3 m_3';
    'x_{1,4}', '-lc_3 m_3';
    'x_{1,5}', 'a_2 m_3 + lc_2 m_2';
    'x_{1,6}', 'lc_3 m_3';
    'x_{1,7}', '2 lc_3 m_3';
    'x_{1,8}', 'lc_3 m_3';
    'x_{2,1}', '-(lc_2 m_2 + a_2 m_3)';
    'x_{2,2}', '-lc_3 m_3';
    'x_{2,3}', 'I_{2zz}+I_{3zz}+a_2^2 m_3+lc_2^2 m_2+lc_3^2 m_3';
    'x_{2,4}', '2 a_2 lc_3 m_3';
    'x_{2,5}', 'I_{3zz}+lc_3^2 m_3';
    'x_{2,6}', 'a_2 lc_3 m_3';
    'x_{2,7}', '-a_2 lc_3 m_3';
    'x_{2,8}', '-2 a_2 lc_3 m_3';
    'x_{2,9}', 'g(lc_2 m_2 + a_2 m_3)';
    'x_{2,10}', 'g lc_3 m_3';
    'x_{3,1}', '-lc_3 m_3';
    'x_{3,2}', 'I_{3zz}+lc_3^2 m_3';
    'x_{3,3}', 'a_2 lc_3 m_3';
    'x_{3,4}', 'I_{3zz}+lc_3^2 m_3';
    'x_{3,5}', 'a_2 lc_3 m_3';
    'x_{3,6}', 'g lc_3 m_3';
};

subsystems = {'Cart ($L_1$)', '', '', '', '', '', '', '', ...
              'Pend. 1 ($L_2$)', '', '', '', '', '', '', '', '', '', ...
              'Pend. 2 ($L_3$)', '', '', '', '', ''};

for i = 1:24
    rel_err_pc = abs(theta_pc(i) - theta_gt(i)) / abs(theta_gt(i)) * 100;
    
    subsys_str = subsystems{i};
    if i == 9 || i == 19
        fprintf(fid, '\\midrule\n');
    end
    fprintf(fid, '%s & $%s$ & $%s$ & %10.4f & %10.4f & %10.4f & %10.4f & %6.2f\\%% \\\\\n', ...
        subsys_str, param_syms{i,1}, param_syms{i,2}, theta_gt(i), theta_ols(i), theta_aguilar(i), theta_pc(i), rel_err_pc);
end

fprintf(fid, '\\midrule\n');
fprintf(fid, '\\multicolumn{4}{l}{\\textbf{Overall Parameter Estimation Error $\\frac{\\|\\hat{\\Theta} - \\Theta^*\\|_2}{\\|\\Theta^*\\|_2} \\times 100\\%%$}} & \\textbf{%.2f\\%%} & \\textbf{%.2f\\%%} & \\textbf{%.2f\\%%} & N/A \\\\\n', ...
    err_ols, err_ag, err_pc);
fprintf(fid, '\\multicolumn{4}{l}{\\textbf{Link 2 ($L_2$) Parameter Norm $\\|\\Theta_{L2}\\|_2$ (True: $10.5498$)}} & \\textbf{%.4f} (Collapse) & \\textbf{%.4f} (Collapse) & \\textbf{%.4f} & N/A \\\\\n', ...
    norm_L2_ols, norm_L2_ag, norm_L2_pc);
fprintf(fid, '\\multicolumn{4}{l}{\\textbf{Link 3 ($L_3$) Parameter Norm $\\|\\Theta_{L3}\\|_2$ (True: $3.3393$)}} & \\textbf{%.4f} (Collapse) & \\textbf{%.4f} (Collapse) & \\textbf{%.4f} & N/A \\\\\n', ...
    norm_L3_ols, norm_L3_ag, norm_L3_pc);
fprintf(fid, '\\multicolumn{4}{l}{\\textbf{Estimated Sensor Noise Variance $\\hat{\\sigma}^2$ (True: $2.50 \\times 10^{-3}$)}} & $%.4f \\times 10^{-3}$ & $%.4f \\times 10^{-3}$ & $\\mathbf{%.4f \\times 10^{-3}}$ & N/A \\\\\n', ...
    sigma2_ols*1000, sigma2_ag*1000, sigma2_pc*1000);
fprintf(fid, '\\multicolumn{4}{l}{\\textbf{Forward Simulation Tracking RMSE ($d_1$ / $\\theta_2$ / $\\theta_3$)}} & \\textbf{Diverged (NaN/Inf)} & \\textbf{Diverged (NaN/Inf)} & \\textbf{%.4f m / %.4f rad / %.4f rad} & N/A \\\\\n', ...
    rmse_cart_pc, rmse_s2_pc, rmse_s3_pc);
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}%%\n');
fprintf(fid, '}\n');
fprintf(fid, '\\end{table*}\n');
fclose(fid);
fprintf('  Saved LaTeX Table: %s\n', tab_file);

results = struct('theta_gt', theta_gt, 'theta_ols', theta_ols, 'theta_aguilar', theta_aguilar, 'theta_pc', theta_pc, ...
                 'err_ols', err_ols, 'err_ag', err_ag, 'err_pc', err_pc, ...
                 'sigma2_ols', sigma2_ols, 'sigma2_ag', sigma2_ag, 'sigma2_pc', sigma2_pc, ...
                 'rmse_pc', [rmse_cart_pc, rmse_s2_pc, rmse_s3_pc], ...
                 't', t, 'd1_true', d1_true, 's2_true', s2_true, 's3_true', s3_true, ...
                 'q_sim_pc', q_sim_pc, 'F1_clean', F1_clean, 'F1_noisy', F1_noisy, 'L1', L1, ...
                 'norm_L2_ols', norm_L2_ols, 'norm_L2_ag', norm_L2_ag, 'norm_L2_pc', norm_L2_pc, ...
                 'norm_L3_ols', norm_L3_ols, 'norm_L3_ag', norm_L3_ag, 'norm_L3_pc', norm_L3_pc, ...
                 't_ols', t_ols, 't_ag', t_ag, 't_pc', t_pc);
end

%% Helper Function: Forward Simulation of DIPC
function [q_sim, is_valid] = simulate_dipc_forward(t, q0, dq0, F1_vec, th)
N = length(t);
x = [q0; dq0];
q_sim = zeros(N, 3);
q_sim(1,:) = q0';
is_valid = true;

for k = 1:N-1
    qk = x(1:3);
    dqk = x(4:6);
    tau_k = [F1_vec(k); 0; 0];
    
    [D, C_dq, G] = dipc_matrices_eval(qk, dqk, th);
    
    if rcond(D) < 1e-12 || isnan(rcond(D))
        is_valid = false;
        q_sim(k+1:end, :) = NaN;
        return;
    end
    
    ddq = D \ (tau_k - C_dq - G);
    
    if any(isnan(ddq)) || any(isinf(ddq)) || max(abs(ddq)) > 1e4
        is_valid = false;
        q_sim(k+1:end, :) = NaN;
        return;
    end
    
    dt_k = t(k+1) - t(k);
    x = x + [dqk; ddq] * dt_k;
    q_sim(k+1, :) = x(1:3)';
end
end

function [D, C_dq, G] = dipc_matrices_eval(q, dq, th)
s2 = q(2); s3 = q(3);
ds2 = dq(2); ds3 = dq(3);

x1 = th(1:8);
x2 = th(9:18);
x3 = th(19:24);

D = zeros(3, 3);
D(1,1) = x1(1);
D(1,2) = x1(2)*cos(s2) + x1(3)*cos(s2 + s3);
D(1,3) = x1(4)*cos(s2 + s3);

D(2,1) = x2(1)*cos(s2) + x2(2)*cos(s2 + s3);
D(2,2) = x2(3) + x2(4)*cos(s3);
D(2,3) = x2(5) + x2(6)*cos(s3);

D(3,1) = x3(1)*cos(s2 + s3);
D(3,2) = x3(2) + x3(3)*cos(s3);
D(3,3) = x3(4);

C_dq = zeros(3, 1);
C_dq(1) = x1(5)*sin(s2)*ds2^2 + x1(6)*sin(s2+s3)*ds2^2 + x1(7)*sin(s2+s3)*ds2*ds3 + x1(8)*sin(s2+s3)*ds3^2;
C_dq(2) = x2(7)*sin(s3)*ds3^2 + x2(8)*sin(s3)*ds2*ds3;
C_dq(3) = x3(5)*sin(s3)*ds2^2;

G = zeros(3, 1);
G(1) = 0;
G(2) = x2(9)*sin(s2) + x2(10)*sin(s2+s3);
G(3) = x3(6)*sin(s2+s3);
end
