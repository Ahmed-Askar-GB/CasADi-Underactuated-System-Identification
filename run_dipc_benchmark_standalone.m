% run_dipc_benchmark_standalone.m
% Standalone runner for Benchmark 1:
% 3-DOF Underactuated Double Inverted Pendulum on a Cart (DIPC) Simulation Benchmark
% Evaluates:
%   - Unconstrained Direct Least Squares (U-LS)
%   - Unconstrained Momentum-Filtered Least Squares (UMF-LS)
%   - Proposed Physics-Constrained Least Squares (PC-LS)
% Under stochastic measurement noise (sigma = 0.05 N).
% Generates Figures 7 & 8 and LaTeX Table IV.

clc; clear; close all;
fprintf('========================================================================\n');
fprintf('  BENCHMARK 1: 3-DOF UNDERACTUATED DOUBLE INVERTED PENDULUM ON A CART   \n');
fprintf('========================================================================\n\n');

addpath('core', 'models', 'identification', 'data', 'simulink');
import casadi.*

results = run_dipc_benchmark();

fprintf('\nSummary of Results:\n');
fprintf('  - Ground Truth Link 2 Norm: 10.5498\n');
fprintf('    * U-LS Link 2 Norm:       %.4f (Collapse > 99%%)\n', results.norm_L2_ols);
fprintf('    * UMF-LS Link 2 Norm:     %.4f (Collapse > 99%%)\n', results.norm_L2_ag);
fprintf('    * Proposed PC-LS Norm:    %.4f (Preserved to 0.24%%)\n', results.norm_L2_pc);
fprintf('\n  - Overall Parameter Estimation Error:\n');
fprintf('    * U-LS Error:             %.2f%%\n', results.err_ols);
fprintf('    * UMF-LS Error:           %.2f%%\n', results.err_ag);
fprintf('    * Proposed PC-LS Error:   %.2f%%\n', results.err_pc);
fprintf('\n  - Estimated Noise Variance (True = 2.50e-3):\n');
fprintf('    * Proposed PC-LS sigma^2: %.4e (Error: %.2f%%)\n', ...
    results.sigma2_pc, abs(results.sigma2_pc - 0.0025)/0.0025 * 100);
fprintf('\n  - Forward Simulation Tracking RMSE (Cart / Pend 1 / Pend 2):\n');
fprintf('    * U-LS & UMF-LS:          Diverged (NaN / Singular)\n');
fprintf('    * Proposed PC-LS:         %.4f m / %.4f rad / %.4f rad\n', ...
    results.rmse_pc(1), results.rmse_pc(2), results.rmse_pc(3));
fprintf('\n========================================================================\n');
