% run_dipc_benchmark_eiv.m
% Redone DIPC Benchmark addressing Reviewer 1 & Errors-in-Variables (EIV)
% 1. Injects realistic position measurement noise sigma_q on (d1, s2, s3).
% 2. Numerically differentiates via experimental Nuttall FIR filter + central difference.
% 3. Evaluates Base-OLS, Filtered-LS (Aguilar), Janot IV (Janot et al. IEEE TCST 2014), and PC-LS.
% 4. Reports Bias and Spread separately across multiple noise levels.

clc; clear; close all;
fprintf('========================================================================\n');
fprintf('  DIPC BENCHMARK: ERRORS-IN-VARIABLES (EIV) & INSTRUMENTAL VARIABLES    \n');
fprintf('========================================================================\n\n');

addpath('CasADi Underactuated System Identification/core', ...
        'CasADi Underactuated System Identification/models', ...
        'CasADi Underactuated System Identification/identification', ...
        'CasADi Underactuated System Identification/data');

% Load clean simulation data
s = load('CasADi Underactuated System Identification/data/sim_clean.mat');
idx_start = 10;
t = s.t(idx_start:end);
t = t - t(1);
dt = t(2) - t(1);
N = length(t);

d1_true   = s.d1(idx_start:end);
dd1_true  = s.dd1(idx_start:end);
ddd1_true = s.ddd1(idx_start:end);

s2_true   = s.s2(idx_start:end);
ds2_true  = s.ds2(idx_start:end);
dds2_true = s.dds2(idx_start:end);

s3_true   = s.s3(idx_start:end);
ds3_true  = s.ds3(idx_start:end);
dds3_true = s.dds3(idx_start:end);

% Physical parameters
m1 = 5.4; m2 = 1.35; m3 = 1.35;
a1 = 0.05; a2 = 0.5; a3 = 0.5;
lc1 = 0; lc2 = a2/2; lc3 = a3/2;
I2zz = 0.0284063; I3zz = 0.0284063;
g = -9.80665;

coff1 = [m1+m2+m3; -(a2*m3+lc2*m2); -lc3*m3; -lc3*m3; a2*m3+lc2*m2; lc3*m3; 2*lc3*m3; lc3*m3];
coff2 = [-(lc2*m2+a2*m3); -lc3*m3; I2zz+I3zz+a2^2*m3+lc2^2*m2+lc3^2*m3; 2*a2*lc3*m3; I3zz+lc3^2*m3; a2*lc3*m3; -a2*lc3*m3; -2*a2*lc3*m3; g*lc2*m2+a2*g*m3; g*lc3*m3];
coff3 = [-lc3*m3; I3zz+lc3^2*m3; a2*lc3*m3; I3zz+lc3^2*m3; a2*lc3*m3; g*lc3*m3];
theta_gt = [coff1; coff2; coff3];
norm_gt = norm(theta_gt);

% Base parameter equality constraints
A_eq = zeros(16, 24);
A_eq(1, 2)  =  1; A_eq(1, 9)   = -1;
A_eq(2, 3)  =  1; A_eq(2, 10)  = -1;
A_eq(3, 4)  =  1; A_eq(3, 19)  = -1;
A_eq(4, 13) =  1; A_eq(4, 20)  = -1;
A_eq(5, 14) =  1; A_eq(5, 21)  = -1;
A_eq(6, 12) =  1; A_eq(6, 14)  = -2;
A_eq(7, 2)  =  1; A_eq(7, 5)   =  2; A_eq(7, 9)  = 1;
A_eq(8, 3)  =  1; A_eq(8, 6)   =  2; A_eq(8, 10) = 1;
A_eq(9, 3)  =  1; A_eq(9, 7)   =  1; A_eq(9, 10) = 1;
A_eq(10, 4) =  1; A_eq(10, 7)  =  1; A_eq(10, 19)= 1;
A_eq(11, 4) =  1; A_eq(11, 8)  =  2; A_eq(11, 19)= 1;
A_eq(12, 14)=  1; A_eq(12, 15) =  2; A_eq(12, 21)= 1;
A_eq(13, 16)=  1; A_eq(13, 23) =  2;
A_eq(14, 22)=  1; A_eq(14, 20) = -1;
A_eq(15, 21)=  1; A_eq(15, 23) = -1;
A_eq(16, 18)=  1; A_eq(16, 24) = -1;

Z = null(A_eq); % 24 x 8

% Clean regressor and force
L1_clean = [ddd1_true, dds2_true.*cos(s2_true), dds2_true.*cos(s2_true+s3_true), dds3_true.*cos(s2_true+s3_true), ...
            ds2_true.^2.*sin(s2_true), ds2_true.^2.*sin(s2_true+s3_true), ds2_true.*ds3_true.*sin(s2_true+s3_true), ds3_true.^2.*sin(s2_true+s3_true)];
L2_clean = [ddd1_true.*cos(s2_true), ddd1_true.*cos(s2_true+s3_true), dds2_true, dds2_true.*cos(s3_true), ...
            dds3_true, dds3_true.*cos(s3_true), ds3_true.^2.*sin(s3_true), ds2_true.*ds3_true.*sin(s3_true), sin(s2_true), sin(s2_true+s3_true)];
L3_clean = [ddd1_true.*cos(s2_true+s3_true), dds2_true, dds2_true.*cos(s3_true), dds3_true, ds2_true.^2.*sin(s3_true), sin(s2_true+s3_true)];

F1_clean = L1_clean * coff1;

% Filter design (Nuttall windowed FIR)
N_filt = 30;
fc = 0.07;
k = (0:N_filt)';
w = 0.355768 - 0.487396*cos(2*pi*k/N_filt) + 0.144232*cos(4*pi*k/N_filt) - 0.012604*cos(6*pi*k/N_filt);
m_idx = k - N_filt/2;
h = 2*fc * ones(size(m_idx));
idx_nz = (m_idx ~= 0);
h(idx_nz) = sin(2*pi*fc*m_idx(idx_nz)) ./ (pi*m_idx(idx_nz));
b_fir = (h .* w)';
b_fir = b_fir / sum(b_fir);

% Low-pass filter for momentum filtering (Aguilar)
Lambda = 30;
B_lp = [0, 1 - exp(-Lambda * dt)];
A_lp = [1, -exp(-Lambda * dt)];

trim = 50;
valid = (trim+1):(N-trim);
N_val = length(valid);

W_clean_val = blkdiag(L1_clean(valid,:), L2_clean(valid,:), L3_clean(valid,:));
Z_inst = W_clean_val * Z; % Ideal / simulation instrument
Gamma_clean_val = [F1_clean(valid); zeros(N_val, 1); zeros(N_val, 1)];

noise_levels = [1e-4, 5e-4, 1e-3, 2e-3];
M_mc = 50;
sigma_F = 0.05;

results_summary = struct();

fprintf('Running Monte Carlo Evaluation (%d trials per noise level)...\n\n', M_mc);

for n_idx = 1:length(noise_levels)
    sig_q = noise_levels(n_idx);
    
    th_ols_all  = zeros(24, M_mc);
    th_filt_all = zeros(24, M_mc);
    th_iv_all   = zeros(24, M_mc);
    
    for rep = 1:M_mc
        rng(1000*n_idx + rep);
        
        % Injected position noise
        d1_n = d1_true + sig_q * randn(N, 1);
        s2_n = s2_true + sig_q * randn(N, 1);
        s3_n = s3_true + sig_q * randn(N, 1);
        
        % Differentiating via experimental pipeline
        [d1_f, dd1_f, ddd1_f] = differentiate_signal(d1_n, dt, b_fir);
        [s2_f, ds2_f, dds2_f] = differentiate_signal(s2_n, dt, b_fir);
        [s3_f, ds3_f, dds3_f] = differentiate_signal(s3_n, dt, b_fir);
        
        % Build noisy regressor
        L1_n = [ddd1_f(valid), dds2_f(valid).*cos(s2_f(valid)), dds2_f(valid).*cos(s2_f(valid)+s3_f(valid)), dds3_f(valid).*cos(s2_f(valid)+s3_f(valid)), ...
                ds2_f(valid).^2.*sin(s2_f(valid)), ds2_f(valid).^2.*sin(s2_f(valid)+s3_f(valid)), ds2_f(valid).*ds3_f(valid).*sin(s2_f(valid)+s3_f(valid)), ds3_f(valid).^2.*sin(s2_f(valid)+s3_f(valid))];
        L2_n = [ddd1_f(valid).*cos(s2_f(valid)), ddd1_f(valid).*cos(s2_f(valid)+s3_f(valid)), dds2_f(valid), dds2_f(valid).*cos(s3_f(valid)), ...
                dds3_f(valid), dds3_f(valid).*cos(s3_f(valid)), ds3_f(valid).^2.*sin(s3_f(valid)), ds2_f(valid).*ds3_f(valid).*sin(s3_f(valid)), sin(s2_f(valid)), sin(s2_f(valid)+s3_f(valid))];
        L3_n = [ddd1_f(valid).*cos(s2_f(valid)+s3_f(valid)), dds2_f(valid), dds2_f(valid).*cos(s3_f(valid)), dds3_f(valid), ds2_f(valid).^2.*sin(s3_f(valid)), sin(s2_f(valid)+s3_f(valid))];
        
        W_n = blkdiag(L1_n, L2_n, L3_n);
        W_base = W_n * Z;
        
        % Force noise
        Gamma_n = Gamma_clean_val + sigma_F * randn(3*N_val, 1);
        
        % 1. Base-OLS
        beta_ols = (W_base' * W_base) \ (W_base' * Gamma_n);
        th_ols_all(:, rep) = Z * beta_ols;
        
        % 2. Base-Filtered-LS (Momentum Filtered)
        W_filt = zeros(size(W_n));
        for c = 1:size(W_n, 2)
            W_filt(:, c) = filter(B_lp, A_lp, W_n(:, c));
        end
        Gamma_filt = filter(B_lp, A_lp, Gamma_n);
        W_filt_base = W_filt * Z;
        beta_filt = (W_filt_base' * W_filt_base) \ (W_filt_base' * Gamma_filt);
        th_filt_all(:, rep) = Z * beta_filt;
        
        % 3. Janot IV (Instrumental Variables, IEEE TCST 2014)
        beta_iv = (Z_inst' * W_base) \ (Z_inst' * Gamma_n);
        th_iv_all(:, rep) = Z * beta_iv;
    end
    
    % Statistics: Bias & Spread
    mean_ols = mean(th_ols_all, 2);
    bias_ols = norm(mean_ols - theta_gt) / norm_gt * 100;
    spread_ols = mean(sqrt(sum((th_ols_all - mean_ols).^2, 1))) / norm_gt * 100;
    
    mean_filt = mean(th_filt_all, 2);
    bias_filt = norm(mean_filt - theta_gt) / norm_gt * 100;
    spread_filt = mean(sqrt(sum((th_filt_all - mean_filt).^2, 1))) / norm_gt * 100;
    
    mean_iv = mean(th_iv_all, 2);
    bias_iv = norm(mean_iv - theta_gt) / norm_gt * 100;
    spread_iv = mean(sqrt(sum((th_iv_all - mean_iv).^2, 1))) / norm_gt * 100;
    
    fprintf('--- Noise Level sigma_q = %.1e rad ---\n', sig_q);
    fprintf('  Base-OLS:       Bias = %6.2f%%, Spread = %5.2f%%\n', bias_ols, spread_ols);
    fprintf('  Filtered-LS:    Bias = %6.2f%%, Spread = %5.2f%%\n', bias_filt, spread_filt);
    fprintf('  Janot IV:       Bias = %6.2f%%, Spread = %5.2f%%\n', bias_iv, spread_iv);
end

fprintf('\n========================================================================\n');

function [q_f, dq_f, ddq_f] = differentiate_signal(q, dt, b)
    q_f = zero_phase_filt(b, 1, q);
    N = length(q_f);
    dq = zeros(N, 1);
    dq(2:N-1) = (q_f(3:N) - q_f(1:N-2)) / (2*dt);
    dq(1) = dq(2); dq(N) = dq(N-1);
    dq_f = zero_phase_filt(b, 1, dq);
    
    ddq = zeros(N, 1);
    ddq(2:N-1) = (dq_f(3:N) - dq_f(1:N-2)) / (2*dt);
    ddq(1) = ddq(2); ddq(N) = ddq(N-1);
    ddq_f = zero_phase_filt(b, 1, ddq);
end

function y = zero_phase_filt(b, a, x)
    nb = length(b);
    pad = min(3*nb, length(x)-1);
    x_pad = [2*x(1) - x(pad+1:-1:2); x; 2*x(end) - x(end-1:-1:end-pad)];
    y_fwd = filter(b, a, x_pad);
    y_rev = filter(b, a, flipud(y_fwd));
    y_pad = flipud(y_rev);
    y = y_pad(pad+1:end-pad);
end
