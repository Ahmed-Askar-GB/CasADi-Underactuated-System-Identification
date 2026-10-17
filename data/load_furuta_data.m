function [diff_data, filter_data] = load_furuta_data(data_file)
    % load_furuta_data: Ingestion & preprocessing of Furuta pendulum data
    % Self-contained: does not depend on Signal Processing Toolbox
    
    if nargin < 1 || isempty(data_file)
        % Search in priority: same directory as this function, ./data, current folder, parent folders
        candidates = {
            fullfile(fileparts(mfilename('fullpath')), 'furutadataexp.mat'), ...
            fullfile('data', 'furutadataexp.mat'), ...
            'furutadataexp.mat', ...
            fullfile('..', 'data', 'furutadataexp.mat'), ...
            fullfile('..', 'Furuta pendulum', 'furutadataexp.mat')
        };
        found = false;
        for c = 1:length(candidates)
            if exist(candidates{c}, 'file')
                data_file = candidates{c};
                found = true;
                break;
            end
        end
        if ~found
            data_file = fullfile('data', 'furutadataexp.mat');
        end
    end
    
    if ~exist(data_file, 'file')
        error('Data file not found: %s', data_file);
    end
    
    data = load(data_file);
    if isfield(data, 'Posicion')
        Pos = data.Posicion;
    else
        error('Variable Posicion not found in %s', data_file);
    end
    
    l = 20000;
    if length(Pos.time) < l
        l = length(Pos.time);
    end
    
    t = Pos.time(1:l);
    q1 = Pos.signals(1,2).values(1:l, 1);
    q2 = Pos.signals(1,3).values(1:l, 1);
    
    KT = 0.1757; % N*m / A
    KI = 2.0;    % A / V
    u1 = Pos.signals(1,1).values(1:l, 1) * KT * KI;
    u2 = zeros(size(u1));
    
    T = mean(diff(t));
    fs = 1/T;
    
    % --- FIR Filter with Nuttall window (identical to FrutunaPen.mlx) ---
    N = 30;
    fc = 0.07; % normalized frequency (relative to fs/2)
    k = (0:N)';
    w = 0.355768 - 0.487396*cos(2*pi*k/N) + 0.144232*cos(4*pi*k/N) - 0.012604*cos(6*pi*k/N);
    m = k - N/2;
    h = 2*fc * ones(size(m));
    idx = (m ~= 0);
    h(idx) = sin(2*pi*fc*m(idx)) ./ (pi*m(idx));
    b = (h .* w)';
    b = b / sum(b);
    
    % Zero-phase filtering with reflection padding
    qf1 = zero_phase_filter(b, 1, q1);
    qf2 = zero_phase_filter(b, 1, q2);
    uf1 = zero_phase_filter(b, 1, u1);
    
    % --- Central difference velocity ---
    dq1 = zeros(l, 1);
    dq2 = zeros(l, 1);
    for i = 2:l-1
        dq1(i) = (qf1(i+1) - qf1(i-1)) / (2*T);
        dq2(i) = (qf2(i+1) - qf2(i-1)) / (2*T);
    end
    dq1(1) = dq1(2); dq1(l) = dq1(l-1);
    dq2(1) = dq2(2); dq2(l) = dq2(l-1);
    
    % Filter velocities for smooth acceleration estimation
    dqf1 = zero_phase_filter(b, 1, dq1);
    dqf2 = zero_phase_filter(b, 1, dq2);
    
    % --- Central difference acceleration ---
    ddq1 = zeros(l, 1);
    ddq2 = zeros(l, 1);
    for i = 2:l-1
        ddq1(i) = (dqf1(i+1) - dqf1(i-1)) / (2*T);
        ddq2(i) = (dqf2(i+1) - dqf2(i-1)) / (2*T);
    end
    ddq1(1) = ddq1(2); ddq1(l) = ddq1(l-1);
    ddq2(1) = ddq2(2); ddq2(l) = ddq2(l-1);
    ddq1 = zero_phase_filter(b, 1, ddq1);
    ddq2 = zero_phase_filter(b, 1, ddq2);
    
    % Trim edge transient artifacts (100 samples)
    trim = 100;
    valid = (trim+1):(l-trim);
    
    diff_data = struct();
    diff_data.t = t(valid);
    diff_data.q = [qf1(valid), qf2(valid)];
    diff_data.dq = [dq1(valid), dq2(valid)];
    diff_data.dqf = [dqf1(valid), dqf2(valid)];
    diff_data.ddq = [ddq1(valid), ddq2(valid)];
    diff_data.tau = [uf1(valid), u2(valid)];
    diff_data.N = length(valid);
    diff_data.T = T;
    diff_data.fs = fs;
    
    % Also save full filtered signals for the Aguilar method
    filter_data = struct();
    filter_data.t = t;
    filter_data.qf1 = qf1;
    filter_data.qf2 = qf2;
    filter_data.dq1 = dq1;
    filter_data.dq2 = dq2;
    filter_data.u1 = u1;
    filter_data.u2 = u2;
    filter_data.T = T;
    filter_data.l = l;
end

function y = zero_phase_filter(b, a, x)
    % Robust zero-phase forward-backward filter
    nb = length(b);
    pad = min(3*nb, length(x)-1);
    % Reflection padding to eliminate boundary transients
    x_pad = [2*x(1) - x(pad+1:-1:2); x; 2*x(end) - x(end-1:-1:end-pad)];
    y_fwd = filter(b, a, x_pad);
    y_rev = filter(b, a, flipud(y_fwd));
    y_pad = flipud(y_rev);
    y = y_pad(pad+1:end-pad);
end
