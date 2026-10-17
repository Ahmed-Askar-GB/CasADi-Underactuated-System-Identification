% run_dipc_simulation.m
% Executes the MATLAB/Simulink Simscape Multibody model of the 3-DOF Underactuated DIPC
% and generates/updates the simulation trajectory dataset (data/Sim1.mat).

clc; clear; close all;
fprintf('========================================================================\n');
fprintf('  DIPC SIMULINK SIMULATION RUNNER (GENERATE SIM1.MAT)                   \n');
fprintf('========================================================================\n');

% Physical parameters conforming to Table I of manuscript
m1 = 5.40;          % Cart mass [kg]
m2 = 1.35;          % Pendulum link 1 mass [kg]
m3 = 1.35;          % Pendulum link 2 mass [kg]
a1 = 0.05;          % Base offset [m]
a2 = 0.50;          % Pendulum 1 length [m]
a3 = 0.50;          % Pendulum 2 length [m]
lc1 = 0.00;         % Cart CoM [m]
lc2 = a2 / 2.0;     % Pendulum 1 CoM [m] (0.25 m)
lc3 = a3 / 2.0;     % Pendulum 2 CoM [m] (0.25 m)
I2zz = 0.0284063;   % Pendulum 1 inertia [kg*m^2]
I3zz = 0.0284063;   % Pendulum 2 inertia [kg*m^2]
g = -9.80665;       % Gravity [m/s^2]

% Multi-frequency excitation parameters
% F1(t) = 80*sin(2*pi*f1*t) + 40*cos(2*pi*f2*t)
f1 = 1.0;           % Excitation frequency 1 [Hz]
f2 = 2.5;           % Excitation frequency 2 [Hz]

model_name = 'Model_under1';
model_file = fullfile(fileparts(mfilename('fullpath')), [model_name, '.slx']);

if ~exist(model_file, 'file')
    error('Simulink model not found: %s', model_file);
end

fprintf('Loading Simulink model: %s...\n', model_file);
load_system(model_file);

% Configure simulation parameters (fixed step 1 ms, 2.0 seconds)
set_param(model_name, 'StopTime', '2.005');
set_param(model_name, 'SolverType', 'Fixed-step');
set_param(model_name, 'Solver', 'ode4');
set_param(model_name, 'FixedStep', '0.001');

fprintf('Running forward simulation (T = 2.0 s, dt = 0.001 s)...\n');
out = sim(model_name);

output_data_file = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'Sim1.mat');
fprintf('Saving simulation output to %s...\n', output_data_file);
save(output_data_file, 'out');

close_system(model_name, 0);
fprintf('Simulation complete. Output successfully saved.\n');
