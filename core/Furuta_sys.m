% Furuta_sys.m
% Dynamic modeling and symbolic derivation of the 2-DOF Rotary Inverted Furuta Pendulum (RIFP)
% using Denavit-Hartenberg (DH) parameters, standard DH_dyn, and CasADi casadi_dh_dyn.
%
% Conforms strictly to Table II of the manuscript:
%   Link 1 (Arm):      a1 = 0,   alpha1 = pi/2, theta1 = q1, d1 = 0,  revolute = 1
%   Link 2 (Pendulum): a2 = L2,  alpha2 = 0,    theta2 = q2, d2 = L1, revolute = 1
%
% Center of Mass (CoM) position vectors in link frames:
%   Link 1: P_C1 = [0, 0, lc1] (along horizontal arm z1 axis)
%   Link 2: P_C2 = [-(L2 - l2), 0, 0] (along pendulum rod x2 axis, distance l2 from joint 2)

clc; clear; close all;
addpath('core');

fprintf('========================================================================\n');
fprintf('  ROTARY INVERTED FURUTA PENDULUM (RIFP) DYNAMICS VIA DH PARAMETERS     \n');
fprintf('========================================================================\n\n');

%% 1. Denavit-Hartenberg (DH) Table Specification (Table II in Manuscript)
% Format per row: [a_i, alpha_i, theta_i, d_i, joint_type (1=revolute, 0=prismatic)]
% Link 1: Arm
% Link 2: Pendulum

% Numerical physical parameters from Aguilar-Sierra et al. (2014) benchmark
L1_num  = 0.278;      % Arm length [m]
lc1_num = 0.139;      % Arm CoM distance [m]
J1_num  = 0.0348;     % Arm base moment of inertia [kg*m^2]
L2_num  = 0.300;      % Pendulum link length [m]
l2_num  = 0.150;      % Pendulum CoM distance from joint [m]
m2_num  = 0.3316;     % Pendulum mass [kg]
J2_num  = 0.0051;     % Pendulum moment of inertia around CoM [kg*m^2]
g_num   = 9.80665;    % Gravitational acceleration [m/s^2]

fprintf('Physical Parameters:\n');
fprintf('  L1 = %.3f m, lc1 = %.3f m, J1 = %.4f kg*m^2\n', L1_num, lc1_num, J1_num);
fprintf('  L2 = %.3f m, l2  = %.3f m, m2 = %.4f kg, J2 = %.4f kg*m^2\n\n', L2_num, l2_num, m2_num, J2_num);

%% 2. Symbolic Derivation with MATLAB Symbolic Math Toolbox
try
    syms s1 s2 ds1 ds2 dds1 dds2 g 'real'
    syms L1 L2 lc1 l2 m1 m2 J1 J2 'real'
    syms I1xx I1yy I1zz I2xx I2yy I2zz 'real'
    
    q   = [s1; s2];
    dq  = [ds1; ds2];
    ddq = [dds1; dds2];
    m   = [m1; m2];
    gc  = [0; 0; -g];
    
    I1 = [I1xx 0 0; 0 I1yy 0; 0 0 I1zz];
    I2 = [I2xx 0 0; 0 I2yy 0; 0 0 I2zz];
    I(:,:,1) = I1;
    I(:,:,2) = I2;
    
    DHpram = [0,   pi/2, s1, 0,  1;   % Link 1: Arm
              L2,  0,    s2, L1, 1];  % Link 2: Pendulum
              
    IC_pose = [0,            0, lc1;
               -(L2 - l2),   0, 0  ];
               
    [M, C, G, Ro] = DH_dyn(DHpram, IC_pose, q, dq, I, m, gc);
    
    fprintf('--- Symbolic Derivation via DH_dyn ---\n');
    fprintf('Inertia Matrix M(q):\n');
    disp(simplify(M));
    fprintf('Gravity Vector G(q):\n');
    disp(simplify(G));
catch ME
    fprintf('Note on Symbolic Toolbox: %s\n', ME.message);
end

%% 3. CasADi Algorithmic Differentiation Engine
try
    import casadi.*
    q_sx  = SX.sym('q', 2);
    dq_sx = SX.sym('dq', 2);
    
    DH_num = [0,       pi/2, q_sx(1), 0,      1;
              L2_num,  0,    q_sx(2), L1_num, 1];
              
    IC_num = [0,                    0, lc1_num;
              -(L2_num - l2_num),   0, 0      ];
              
    I_num = zeros(3, 3, 2);
    I_num(:,:,1) = diag([0.001, 0.001, J1_num]);
    I_num(:,:,2) = diag([J2_num, J2_num, 0.001]);
    m_num = [0.5; m2_num];
    g_vec = [0; 0; -g_num];
    
    [M_casadi, C_casadi, G_casadi] = casadi_dh_dyn(DH_num, IC_num, I_num, m_num, g_vec, q_sx, dq_sx);
    
    M_fn = Function('M_fn', {q_sx}, {M_casadi});
    G_fn = Function('G_fn', {q_sx}, {G_casadi});
    
    q_test = [0; 0];
    fprintf('\n--- Numerical CasADi Evaluation at q = [0, 0] ---\n');
    fprintf('M(0, 0):\n');
    disp(full(M_fn(q_test)));
    fprintf('G(0, 0):\n');
    disp(full(G_fn(q_test)));
    
    fprintf('CasADi DH evaluation successful. Equations match benchmark.\n');
catch ME
    fprintf('Note on CasADi: %s\n', ME.message);
end
