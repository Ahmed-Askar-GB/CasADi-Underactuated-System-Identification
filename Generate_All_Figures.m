% Generate_All_Figures.m
% Standalone Publication Figure Generator
% Framework: Physics-Constrained Data-Driven Modeling via CasADi for Underactuated Systems
% 
% PURPOSE:
%   This script loads the pre-computed identification and simulation results
%   from 'data/identification_results.mat' and generates all publication-quality
%   figures (both vector .pdf and 300-DPI .png) without requiring re-running
%   any optimization or identification solvers.
%
% USAGE:
%   Simply run `Generate_All_Figures` in MATLAB.
%   To customize appearance (colors, font sizes, line widths, export formats),
%   modify the parameters in the "STYLING & CONFIGURATION" section below.

clc; clear; close all;

fprintf('========================================================================\n');
fprintf('  PUBLICATION FIGURE GENERATOR (WITHOUT RE-RUNNING OPTIMIZATION)        \n');
fprintf('========================================================================\n');

%% ========================================================================
%% 1. STYLING & CONFIGURATION OPTIONS (USER EDITABLE)
%% ========================================================================

% Figures to generate (select any subset from 1 to 8, or 1:8 for all)
%   1: RIFP Experimental Trajectory Overview
%   2: RIFP Torque Tracking & Reconstruction Residuals
%   3: RIFP Parameter Recovery Ratio & Shrinking Window Robustness
%   4: RIFP Multi-Step Dynamic Forward Simulation Validation
%   5: RIFP Inertia Positive Definiteness min det(D) vs Time
%   6: RIFP Nonlinear Regularized Friction Model Characteristic Curves
%   7: DIPC Subsystem Relative Parameter Error Bar Chart
%   8: DIPC Forward Simulation Tracking Verification
figs_to_plot = 1:8;

% Output Directories & Formats
if exist('paper', 'dir')
    fig_dir = fullfile('paper', 'figures');
else
    fig_dir = 'figures';
end
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

export_pdf      = true;     % Export vector PDF
export_png      = true;     % Export high-res PNG
png_resolution  = 300;      % DPI for PNG exports

% Typography & Global Font Settings
font_name           = 'Times New Roman';
font_size_axes      = 8.5;
font_size_labels    = 9.0;
font_size_titles    = 9.5;
font_size_legends   = 7.5;

% Line Widths & Marker Sizes
lw_bold     = 1.4;
lw_main     = 1.2;
lw_sub      = 1.0;
lw_thin     = 0.8;
ms_main     = 5.0;

% Color Palette (Professional Publication Standards)
col_gt          = [0.15, 0.15, 0.15];   % Ground Truth / Measured: Jet Black / Charcoal
col_prop        = [0.00, 0.45, 0.74];   % Proposed FBC-LS (or PC-LS): Deep Blue
col_baseline    = [0.85, 0.10, 0.15];   % UMF-LS Baseline (Aguilar): Crimson / Red
col_direct      = [0.18, 0.55, 0.34];   % U-LS: Forest Green
col_decoupled   = [0.90, 0.32, 0.00];   % D-LS: Deep Orange
col_gray_pts    = [0.65, 0.65, 0.65];   % Noisy Scatter: Light Gray

%% ========================================================================
%% 2. LOAD PRE-COMPUTED BENCHMARK DATA
%% ========================================================================

candidates = {
    fullfile('data', 'identification_results.mat'), ...
    fullfile(fileparts(mfilename('fullpath')), 'data', 'identification_results.mat'), ...
    'identification_results.mat'
};
data_file = '';
for c = 1:length(candidates)
    if exist(candidates{c}, 'file')
        data_file = candidates{c};
        break;
    end
end

if isempty(data_file)
    error(['Data file "data/identification_results.mat" not found.\n', ...
           'Please execute "main_identification" once to generate the workspace results.']);
end

fprintf('[Step 1/2] Loading pre-computed workspace from %s...\n', data_file);
loaded_data = load(data_file);
fprintf('  Data loaded successfully.\n');

% Apply Global Figure Properties
set(groot, 'DefaultAxesFontName', font_name);
set(groot, 'DefaultTextFontName', font_name);
set(groot, 'DefaultAxesFontSize', font_size_axes);
set(groot, 'DefaultLineLineWidth', lw_main);

fprintf('[Step 2/2] Generating Publication Figures...\n\n');

%% ========================================================================
%% FIGURE 1: RIFP Experimental Trajectory Overview
%% ========================================================================
if ismember(1, figs_to_plot)
    fprintf('  -> Generating Figure 1: Experimental Trajectories (RIFP)...\n');
    f1 = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 6.8, 3.6], 'Color', 'w');
    t1 = tiledlayout(f1, 3, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
    
    t_data = loaded_data.t;
    q_data = loaded_data.q;
    tau_data = loaded_data.tau;
    
    % Tile 1: Arm Angle q1
    nexttile;
    plot(t_data, q_data(:,1), 'Color', col_prop, 'LineWidth', lw_main); hold on; grid on; box on;
    set(gca, 'FontSize', font_size_axes);
    ylabel('Arm $q_1$ [rad]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    xlim([t_data(1), t_data(end)]); xticklabels([]);
    ylim([-3.5, 3.5]);
    
    % Tile 2: Pendulum Angle q2
    nexttile;
    plot(t_data, q_data(:,2), 'Color', col_baseline, 'LineWidth', lw_main); hold on; grid on; box on;
    set(gca, 'FontSize', font_size_axes);
    ylabel('Pendulum $q_2$ [rad]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    xlim([t_data(1), t_data(end)]); xticklabels([]);
    ylim([-3.5, 3.5]);
    
    % Tile 3: Control Input Torque tau1
    nexttile;
    plot(t_data, tau_data(:,1), 'Color', col_gt, 'LineWidth', lw_sub); hold on; grid on; box on;
    set(gca, 'FontSize', font_size_axes);
    ylabel('Torque $\tau_1$ [N$\cdot$m]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    xlim([t_data(1), t_data(end)]);
    ylim([-0.35, 0.35]);
    
    if export_pdf
        exportgraphics(f1, fullfile(fig_dir, 'fig1_experimental_trajectories.pdf'), 'ContentType', 'vector');
    end
    if export_png
        exportgraphics(f1, fullfile(fig_dir, 'fig1_experimental_trajectories.png'), 'Resolution', png_resolution);
    end
    close(f1);
    fprintf('     Saved fig1_experimental_trajectories (.pdf / .png)\n');
end

%% ========================================================================
%% FIGURE 2: RIFP Dynamic Torque Reconstruction Comparison
%% ========================================================================
if ismember(2, figs_to_plot)
    fprintf('  -> Generating Figure 2: Torque Tracking & Residuals (RIFP)...\n');
    f2 = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 6.8, 3.2], 'Color', 'w');
    t2 = tiledlayout(f2, 2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
    
    t_zoom = loaded_data.t_zoom;
    tau_meas_zoom = loaded_data.tau_meas_zoom;
    tau_pcls_zoom = loaded_data.tau_pcls_zoom;
    tau_ag_zoom   = loaded_data.tau_ag_zoom;
    tau_ls_zoom   = loaded_data.tau_ls_zoom;
    
    % Tile 1: Torque Tracking
    nexttile;
    plot(t_zoom, tau_meas_zoom, '-', 'Color', col_gt, 'LineWidth', lw_bold, 'DisplayName', 'Measured'); hold on;
    plot(t_zoom, tau_pcls_zoom, '--', 'Color', col_prop, 'LineWidth', lw_main, 'DisplayName', 'FBC-LS (Proposed)');
    plot(t_zoom, tau_ag_zoom, ':', 'Color', col_baseline, 'LineWidth', lw_main, 'DisplayName', 'UMF-LS Baseline');
    plot(t_zoom, tau_ls_zoom, '-.', 'Color', col_direct, 'LineWidth', lw_thin, 'DisplayName', 'U-LS');
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    ylabel('Torque [N$\cdot$m]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'Orientation', 'horizontal', 'FontSize', font_size_legends);
    xlim([t_zoom(1), t_zoom(end)]); xticklabels([]);
    ylim([-0.32, 0.35]);
    
    % Tile 2: Reconstruction Residuals
    nexttile;
    plot(t_zoom, tau_meas_zoom - tau_pcls_zoom, '-', 'Color', col_prop, 'LineWidth', lw_sub, 'DisplayName', 'FBC-LS Residual'); hold on;
    plot(t_zoom, tau_meas_zoom - tau_ag_zoom, ':', 'Color', col_baseline, 'LineWidth', lw_sub, 'DisplayName', 'UMF-LS Residual');
    plot(t_zoom, tau_meas_zoom - tau_ls_zoom, '-.', 'Color', col_direct, 'LineWidth', lw_thin, 'DisplayName', 'U-LS Residual');
    yline(0, 'k-', 'LineWidth', lw_thin);
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    ylabel('Residual [N$\cdot$m]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    legend('Location', 'southeast', 'Orientation', 'horizontal', 'FontSize', font_size_legends);
    xlim([t_zoom(1), t_zoom(end)]);
    ylim([-0.22, 0.22]);
    
    if export_pdf
        exportgraphics(f2, fullfile(fig_dir, 'fig2_torque_tracking_comparison.pdf'), 'ContentType', 'vector');
    end
    if export_png
        exportgraphics(f2, fullfile(fig_dir, 'fig2_torque_tracking_comparison.png'), 'Resolution', png_resolution);
    end
    close(f2);
    fprintf('     Saved fig2_torque_tracking_comparison (.pdf / .png)\n');
end

%% ========================================================================
%% FIGURE 3: Parameter Recovery & Robustness Demonstration (2-Panel)
%% ========================================================================
if ismember(3, figs_to_plot)
    fprintf('  -> Generating Figure 3: Parameter Recovery & Window Robustness...\n');
    f3 = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 7.0, 4.7], 'Color', 'w');
    t3 = tiledlayout(f3, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Panel (a): Normalized Parameter Recovery Ratio
    nexttile;
    th_dec  = loaded_data.th_decoupled;
    th_dls  = loaded_data.th_direct_ls;
    th_ag   = loaded_data.th_aguilar;
    th_pc   = loaded_data.th_pcls_prop;
    
    ratio_data = [th_dec ./ th_pc, th_dls ./ th_pc, th_ag ./ th_pc, th_pc ./ th_pc];
    b_bar = bar(1:9, ratio_data, 0.85);
    b_bar(1).FaceColor = col_decoupled;
    b_bar(2).FaceColor = [0.98, 0.75, 0.18]; % Gold/Yellow (U-LS)
    b_bar(3).FaceColor = col_direct;          % Forest Green (UMF-LS)
    b_bar(4).FaceColor = col_prop;            % Deep Blue (FBC-LS)
    yline(1.0, 'k--', 'LineWidth', 1.0);
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    set(gca, 'XTick', 1:9, 'XTickLabel', ...
        {'$\theta_1$', '$\theta_2$', '$\theta_3$', '$\theta_4$', '$\theta_5$', '$f_{v1}$', '$f_{v2}$', '$f_{c1}$', '$f_{c2}$'}, ...
        'TickLabelInterpreter', 'latex', 'FontSize', font_size_axes);
    ylabel('Normalized $\hat{\theta}/\theta_{\mathrm{ref}}$', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    ylim([0, 2.65]);
    title('(a) Parameter Recovery Across Methodologies Normalized to Physical Reference', 'FontWeight', 'bold', 'FontSize', font_size_titles);
    legend({'D-LS (Null Collapse)', 'U-LS', 'UMF-LS Baseline', 'FBC-LS (Proposed)'}, ...
        'Location', 'northeast', 'Orientation', 'horizontal', 'FontSize', font_size_legends);
    
    % Panel (b): Observation Window Duration Robustness (Table VI data)
    nexttile;
    T_win = [0.3, 0.5, 0.8, 1.0, 2.0];
    det_mfls_x1e4 = [-0.117, -0.367, 1.68, 1.85, 4.06];
    det_pcls_x1e4 = [1.00, 1.00, 2.93, 3.45, 4.06];
    p1 = plot(T_win, det_mfls_x1e4, 's--', 'Color', col_decoupled, 'LineWidth', lw_bold, 'MarkerSize', ms_main); hold on;
    p2 = plot(T_win, det_pcls_x1e4, 'o-', 'Color', col_prop, 'LineWidth', lw_bold, 'MarkerSize', ms_main);
    yline(0, 'k-', 'LineWidth', lw_thin);
    yline(1.0, ':', 'Color', col_prop, 'LineWidth', lw_thin);
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    xlim([0.2, 2.1]);
    ylim([-0.6, 4.8]);
    xlabel('Observation Horizon Duration $T$ [s]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    ylabel('$\min \det(D) \times 10^4$', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold', 'Color', col_prop);
    title('(b) Physical Validity and Forward Simulation Stability Under Shrinking Observation Horizons', 'FontWeight', 'bold', 'FontSize', font_size_titles);
    
    % Twin y-axis for forward simulation RMSE
    yyaxis right;
    rmse_mfls_data = [NaN, 0.9377, 0.2140, 0.1852, 0.0152];
    rmse_pcls_data = [0.0784, 0.0617, 0.0381, 0.0246, 0.0114];
    p3 = plot(T_win, rmse_mfls_data, '^--', 'Color', col_baseline, 'LineWidth', lw_main, 'MarkerSize', ms_main);
    p4 = plot(T_win, rmse_pcls_data, 'd-', 'Color', col_direct, 'LineWidth', lw_main, 'MarkerSize', ms_main);
    ylabel('Sim RMSE $q_2$ [rad]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold', 'Color', col_baseline);
    ylim([0, 1.1]);
    legend([p1, p2, p3, p4], {'UMF-LS $\min\det(D)\times 10^4$', 'FBC-LS $\min\det(D)\times 10^4$', 'UMF-LS Sim RMSE $q_2$', 'FBC-LS Sim RMSE $q_2$'}, ...
        'Location', 'southeast', 'FontSize', font_size_legends, 'NumColumns', 2);
    
    if export_pdf
        exportgraphics(f3, fullfile(fig_dir, 'fig3_nullspace_collapse_demonstration.pdf'), 'ContentType', 'vector');
    end
    if export_png
        exportgraphics(f3, fullfile(fig_dir, 'fig3_nullspace_collapse_demonstration.png'), 'Resolution', png_resolution);
    end
    close(f3);
    fprintf('     Saved fig3_nullspace_collapse_demonstration (.pdf / .png)\n');
end

%% ========================================================================
%% FIGURE 4: RIFP Multi-Step Dynamic Forward Simulation Validation
%% ========================================================================
if ismember(4, figs_to_plot)
    fprintf('  -> Generating Figure 4: Forward Simulation Validation (RIFP)...\n');
    f4 = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 6.8, 3.2], 'Color', 'w');
    t4 = tiledlayout(f4, 2, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
    
    t_val = loaded_data.t_val;
    q_val = loaded_data.q_val;
    q_sim_pcls = loaded_data.q_sim_pcls;
    q_sim_ag   = loaded_data.q_sim_ag;
    q_sim_dls  = loaded_data.q_sim_dls;
    
    % Tile 1: Arm Angle q1
    nexttile;
    plot(t_val, q_val(:,1), '-', 'Color', col_gt, 'LineWidth', lw_bold, 'DisplayName', 'Ground Truth'); hold on;
    plot(t_val, q_sim_pcls(:,1), '--', 'Color', col_prop, 'LineWidth', lw_main, 'DisplayName', 'FBC-LS (Proposed)');
    plot(t_val, q_sim_ag(:,1), ':', 'Color', col_baseline, 'LineWidth', lw_main, 'DisplayName', 'UMF-LS Baseline');
    plot(t_val, q_sim_dls(:,1), '-.', 'Color', col_direct, 'LineWidth', lw_thin, 'DisplayName', 'U-LS');
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    ylabel('Arm $q_1$ [rad]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'Orientation', 'horizontal', 'FontSize', font_size_legends);
    xlim([t_val(1), t_val(end)]); xticklabels([]);
    ylim([-2.0, 3.2]);
    
    % Tile 2: Pendulum Angle q2
    nexttile;
    plot(t_val, q_val(:,2), '-', 'Color', col_gt, 'LineWidth', lw_bold, 'DisplayName', 'Ground Truth'); hold on;
    plot(t_val, q_sim_pcls(:,2), '--', 'Color', col_prop, 'LineWidth', lw_main, 'DisplayName', 'FBC-LS (Proposed)');
    plot(t_val, q_sim_ag(:,2), ':', 'Color', col_baseline, 'LineWidth', lw_main, 'DisplayName', 'UMF-LS Baseline');
    plot(t_val, q_sim_dls(:,2), '-.', 'Color', col_direct, 'LineWidth', lw_thin, 'DisplayName', 'U-LS');
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    ylabel('Pendulum $q_2$ [rad]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    xlim([t_val(1), t_val(end)]);
    ylim([-3.2, 3.2]);
    
    if export_pdf
        exportgraphics(f4, fullfile(fig_dir, 'fig4_forward_simulation_validation.pdf'), 'ContentType', 'vector');
    end
    if export_png
        exportgraphics(f4, fullfile(fig_dir, 'fig4_forward_simulation_validation.png'), 'Resolution', png_resolution);
    end
    close(f4);
    fprintf('     Saved fig4_forward_simulation_validation (.pdf / .png)\n');
end

%% ========================================================================
%% FIGURE 5: RIFP Inertia Positive Definiteness min det(D) vs Time
%% ========================================================================
if ismember(5, figs_to_plot)
    fprintf('  -> Generating Figure 5: Inertia Positive Definiteness...\n');
    f5 = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 6.0, 2.5], 'Color', 'w');
    
    t_data = loaded_data.t;
    q_data = loaded_data.q;
    th_pc  = loaded_data.th_pcls_prop;
    th_ag  = loaded_data.th_aguilar;
    th_dls = loaded_data.th_direct_ls;
    
    det_traj_pcls = (th_pc(1) + th_pc(2)*(sin(q_data(:,2)).^2)) * th_pc(4) - (th_pc(3)^2)*(cos(q_data(:,2)).^2);
    det_traj_ag   = (th_ag(1) + th_ag(2)*(sin(q_data(:,2)).^2)) * th_ag(4) - (th_ag(3)^2)*(cos(q_data(:,2)).^2);
    det_traj_dls  = (th_dls(1) + th_dls(2)*(sin(q_data(:,2)).^2)) * th_dls(4) - (th_dls(3)^2)*(cos(q_data(:,2)).^2);
    
    plot(t_data, det_traj_pcls * 1e4, '-', 'Color', col_prop, 'LineWidth', lw_bold, 'DisplayName', 'FBC-LS $\det(D) \times 10^4$'); hold on;
    plot(t_data, det_traj_ag * 1e4, '--', 'Color', col_baseline, 'LineWidth', lw_main, 'DisplayName', 'UMF-LS $\det(D) \times 10^4$');
    plot(t_data, det_traj_dls * 1e4, '-.', 'Color', col_direct, 'LineWidth', lw_sub, 'DisplayName', 'U-LS $\det(D) \times 10^4$');
    yline(0, 'k--', 'LineWidth', lw_thin, 'DisplayName', 'Singularity Boundary ($\det=0$)');
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    ylabel('$\det(D(\mathbf{q})) \times 10^4$ [$\mathrm{kg}^2\cdot\mathrm{m}^4$]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', font_size_legends);
    xlim([t_data(1), t_data(end)]);
    ylim([-0.5, 8.5]);
    
    if export_pdf
        exportgraphics(f5, fullfile(fig_dir, 'fig5_physics_constraints_validation.pdf'), 'ContentType', 'vector');
    end
    if export_png
        exportgraphics(f5, fullfile(fig_dir, 'fig5_physics_constraints_validation.png'), 'Resolution', png_resolution);
    end
    close(f5);
    fprintf('     Saved fig5_physics_constraints_validation (.pdf / .png)\n');
end

%% ========================================================================
%% FIGURE 6: RIFP Nonlinear Regularized Friction Model Curves
%% ========================================================================
if ismember(6, figs_to_plot)
    fprintf('  -> Generating Figure 6: Nonlinear Friction Characteristic Curves...\n');
    f6 = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 6.5, 2.5], 'Color', 'w');
    t6 = tiledlayout(f6, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    v_eval  = loaded_data.v_eval;
    F1_sym  = loaded_data.F1_sym;
    F1_asym = loaded_data.F1_asym;
    F2_sym  = loaded_data.F2_sym;
    
    % Tile 1: Arm Friction
    nexttile;
    plot(v_eval, F1_sym, '-', 'Color', col_prop, 'LineWidth', lw_main, 'DisplayName', 'Symmetric Model'); hold on;
    plot(v_eval, F1_asym, '--', 'Color', col_baseline, 'LineWidth', lw_main, 'DisplayName', 'Asymmetric Model');
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    xlabel('Arm Velocity $\dot{q}_1$ [rad/s]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    ylabel('Friction Torque $F_1$ [N$\cdot$m]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    legend('Location', 'northwest', 'FontSize', font_size_legends);
    xlim([-4, 4]);
    
    % Tile 2: Pendulum Friction
    nexttile;
    plot(v_eval, F2_sym, '-', 'Color', col_prop, 'LineWidth', lw_main, 'DisplayName', 'Symmetric Model'); hold on;
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    xlabel('Pendulum Velocity $\dot{q}_2$ [rad/s]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    ylabel('Friction Torque $F_2$ [N$\cdot$m]', 'Interpreter', 'latex', 'FontSize', font_size_labels, 'FontWeight', 'bold');
    legend('Location', 'northwest', 'FontSize', font_size_legends);
    xlim([-4, 4]);
    
    if export_pdf
        exportgraphics(f6, fullfile(fig_dir, 'fig6_friction_characteristic_curves.pdf'), 'ContentType', 'vector');
    end
    if export_png
        exportgraphics(f6, fullfile(fig_dir, 'fig6_friction_characteristic_curves.png'), 'Resolution', png_resolution);
    end
    close(f6);
    fprintf('     Saved fig6_friction_characteristic_curves (.pdf / .png)\n');
end

%% ========================================================================
%% FIGURE 7 (DIPC Figure A): Relative Parameter Estimation Error
%% ========================================================================
if ismember(7, figs_to_plot) && isfield(loaded_data, 'dipc_results')
    fprintf('  -> Generating Figure 7: Subsystem Relative Parameter Error (DIPC)...\n');
    dipc = loaded_data.dipc_results;
    
    err_bar_ols = abs(dipc.theta_ols - dipc.theta_gt) ./ abs(dipc.theta_gt) * 100;
    err_bar_ag  = abs(dipc.theta_aguilar - dipc.theta_gt) ./ abs(dipc.theta_gt) * 100;
    err_bar_pc  = abs(dipc.theta_pc - dipc.theta_gt) ./ abs(dipc.theta_gt) * 100;
    bar_data = [err_bar_ols, err_bar_ag, err_bar_pc];
    
    figA = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1, 1, 6.8, 4.4]);
    tA = tiledlayout(figA, 3, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
    
    % Tile 1: Link 1 (Cart, 8 parameters)
    nexttile;
    b1 = bar(1:8, bar_data(1:8, :), 0.85);
    b1(1).FaceColor = [0.85, 0.35, 0.35]; b1(1).EdgeColor = 'none'; % Crimson
    b1(2).FaceColor = col_direct;          b1(2).EdgeColor = 'none'; % Green
    b1(3).FaceColor = col_prop;            b1(3).EdgeColor = 'none'; % Deep Blue
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    xlim([0.4, 8.6]);
    ylim([0, 0.52]);
    xticks(1:8);
    xticklabels({'x_{1,1}', 'x_{1,2}', 'x_{1,3}', 'x_{1,4}', 'x_{1,5}', 'x_{1,6}', 'x_{1,7}', 'x_{1,8}'});
    ylabel({'Cart ($L_1$)', 'Rel. Err. [\%]'}, 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', font_size_axes);
    legend({'U-LS', 'UMF-LS Baseline', 'PC-LS (Proposed)'}, ...
        'Location', 'northeast', 'Orientation', 'horizontal', 'FontSize', font_size_legends);
    
    % Tile 2: Link 2 (Passive Pendulum 1, 10 parameters)
    nexttile;
    b2 = bar(1:10, bar_data(9:18, :), 0.85);
    b2(1).FaceColor = [0.85, 0.35, 0.35]; b2(1).EdgeColor = 'none';
    b2(2).FaceColor = col_direct;          b2(2).EdgeColor = 'none';
    b2(3).FaceColor = col_prop;            b2(3).EdgeColor = 'none';
    yline(100, 'k:', 'LineWidth', 0.8);
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    xlim([0.4, 10.6]);
    ylim([0, 115]);
    xticks(1:10);
    xticklabels({'x_{2,1}', 'x_{2,2}', 'x_{2,3}', 'x_{2,4}', 'x_{2,5}', 'x_{2,6}', 'x_{2,7}', 'x_{2,8}', 'x_{2,9}', 'x_{2,10}'});
    ylabel({'Pend. 1 ($L_2$)', 'Rel. Err. [\%]'}, 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', font_size_axes);
    
    % Tile 3: Link 3 (Passive Pendulum 2, 6 parameters)
    nexttile;
    b3 = bar(1:6, bar_data(19:24, :), 0.85);
    b3(1).FaceColor = [0.85, 0.35, 0.35]; b3(1).EdgeColor = 'none';
    b3(2).FaceColor = col_direct;          b3(2).EdgeColor = 'none';
    b3(3).FaceColor = col_prop;            b3(3).EdgeColor = 'none';
    yline(100, 'k:', 'LineWidth', 0.8);
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    xlim([0.4, 6.6]);
    ylim([0, 115]);
    xticks(1:6);
    xticklabels({'x_{3,1}', 'x_{3,2}', 'x_{3,3}', 'x_{3,4}', 'x_{3,5}', 'x_{3,6}'});
    ylabel({'Pend. 2 ($L_3$)', 'Rel. Err. [\%]'}, 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', font_size_axes);
    xlabel('Identified Dynamic Base Parameters', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', font_size_labels);
    
    if export_pdf
        exportgraphics(figA, fullfile(fig_dir, 'fig_dipc_nullspace.pdf'), 'ContentType', 'vector');
    end
    if export_png
        exportgraphics(figA, fullfile(fig_dir, 'fig_dipc_nullspace.png'), 'Resolution', png_resolution);
    end
    close(figA);
    fprintf('     Saved fig_dipc_nullspace (.pdf / .png)\n');
end

%% ========================================================================
%% FIGURE 8 (DIPC Figure B): Forward Simulation Tracking Verification
%% ========================================================================
if ismember(8, figs_to_plot) && isfield(loaded_data, 'dipc_results')
    fprintf('  -> Generating Figure 8: Forward Simulation Tracking (DIPC)...\n');
    dipc = loaded_data.dipc_results;
    
    t_dipc = dipc.t;
    d1_true = dipc.d1_true;
    s2_true = dipc.s2_true;
    s3_true = dipc.s3_true;
    q_sim_pc = dipc.q_sim_pc;
    F1_clean = dipc.F1_clean;
    F1_noisy = dipc.F1_noisy;
    L1 = dipc.L1;
    theta_pc = dipc.theta_pc;
    
    figB = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1, 1, 6.8, 4.4]);
    tB = tiledlayout(figB, 4, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
    
    % Tile 1: Cart Position d1
    nexttile;
    plot(t_dipc, d1_true, 'Color', col_gt, 'LineWidth', lw_bold); hold on;
    plot(t_dipc, q_sim_pc(:,1), '--', 'Color', col_prop, 'LineWidth', lw_main);
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    ylabel('Cart $d_1$ [m]', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', font_size_labels);
    xlim([t_dipc(1), t_dipc(end)]); xticklabels([]);
    ylim([-0.02, 0.40]);
    legend({'Ground Truth', 'PC-LS (Proposed)'}, 'Location', 'northwest', 'FontSize', font_size_legends);
    
    % Tile 2: First Pendulum Angle theta2
    nexttile;
    plot(t_dipc, s2_true, 'Color', col_gt, 'LineWidth', lw_bold); hold on;
    plot(t_dipc, q_sim_pc(:,2), '--', 'Color', col_prop, 'LineWidth', lw_main);
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    ylabel('$\theta_2$ [rad]', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', font_size_labels);
    xlim([t_dipc(1), t_dipc(end)]); xticklabels([]);
    
    % Tile 3: Second Pendulum Angle theta3
    nexttile;
    plot(t_dipc, s3_true, 'Color', col_gt, 'LineWidth', lw_bold); hold on;
    plot(t_dipc, q_sim_pc(:,3), '--', 'Color', col_prop, 'LineWidth', lw_main);
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    ylabel('$\theta_3$ [rad]', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', font_size_labels);
    xlim([t_dipc(1), t_dipc(end)]); xticklabels([]);
    
    % Tile 4: Cart Actuation Force F1
    nexttile;
    plot(t_dipc, F1_clean, 'Color', col_gt, 'LineWidth', lw_main); hold on;
    plot(t_dipc, L1 * theta_pc(1:8), '--', 'Color', col_prop, 'LineWidth', lw_sub);
    grid on; box on; set(gca, 'FontSize', font_size_axes);
    xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', font_size_labels);
    ylabel('Force $F_1$ [N]', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', font_size_labels);
    xlim([t_dipc(1), t_dipc(end)]);
    ylim([-95, 185]);
    
    if export_pdf
        exportgraphics(figB, fullfile(fig_dir, 'fig_dipc_tracking.pdf'), 'ContentType', 'vector');
    end
    if export_png
        exportgraphics(figB, fullfile(fig_dir, 'fig_dipc_tracking.png'), 'Resolution', png_resolution);
    end
    close(figB);
    fprintf('     Saved fig_dipc_tracking (.pdf / .png)\n');
end

fprintf('\n========================================================================\n');
fprintf('  ALL REQUESTED FIGURES GENERATED SUCCESSFULLY IN %s\n', fig_dir);
fprintf('========================================================================\n');
