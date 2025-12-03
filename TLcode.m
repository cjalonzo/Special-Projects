% Transmission-line ladder time-domain simulator (lumped N-section)


clear; close all; clc;

% ---- Physical / line parameters ----
len = 1.0;            % total line length (m)
N = 200;              % number of sections (higher => closer to continuous)
dx = len / N;

% per-unit-length parameters (example values; adjust as needed)
R_per_m = 1.0;        % series resistance (ohm/m)
L_per_m = 250e-9;     % series inductance (H/m)
G_per_m = 1e-3;       % shunt conductance (S/m)
C_per_m = 100e-12;    % shunt capacitance (F/m)


% section (lumped) values
R = R_per_m * dx;     % series resistance per section
L = L_per_m * dx;     % series inductance per section
G = G_per_m * dx;     % shunt conductance per section
C = C_per_m * dx;     % shunt capacitance per section

% --- Skin effect modeling coefficient (adjust as needed) ---
% R_skin scales how strongly resistance increases at high freq.
R_skin_coeff = 0.7;      % empirical factor (0.3–1.5 typical)
f_equiv = 1/(4*dt);      % approximate effective freq content
R_skin = R_skin_coeff * sqrt(f_equiv);

% Effective series resistance including skin effect
R_eff = R + R_skin;


% ---- simulation parameters ----
fs = 5e9;             % sampling frequency (Hz) - choose high enough
dt = 1/fs;
tmax = 200e-9;
steps = floor(tmax / dt);

% Source & load
Vs = 1.0;             % step/pulse amplitude (V)
Rs = 50;              % source series resistance (ohm)
RL = 75;              % load resistance at far end (ohm)

% Preallocate
v = zeros(N,1);       % capacitor voltages at nodes 1..N
i = zeros(N+1,1);     % series currents between nodes 0..N (i(1) = between source and node1)
tvec = (0:steps-1)*dt;

% Gaussian pulse parameters
t0 = 30e-9;
tau = 8e-9;
pulse = @(t) Vs * exp(-((t - t0)/tau).^2);

% Record nodes to observe
obs_idx = [1, round(N/4), round(N/2), N]; % positions: near source, quarter, half, load
Vrec = zeros(length(obs_idx), steps);

% Precompute helper constants (explicit update scheme)
% Use Euler forward for updates (stable if dt small)
for n = 1:steps
    t = (n-1)*dt;
    % update currents (i between nodes): i_k^{n+1} = i_k^{n} + (dt/L)*(v_{k-1} - v_k - R*i_k)
    % note: v_0 is source node voltage (after Rs)
    v_source = (pulse(t) - Rs * i(1));  % approximate source node voltage (voltage division)
    
    % update all internal currents
    for k = 1:N
        v_left = (k==1) * v_source + (k>1) * v(k-1);
        v_right = v(k);
        di = (dt / L) * (v_left - v_right - R_eff * i(k+1)); 
        % note: i array offset chosen so i(k+1) corresponds to series current between node k and k+1
        i(k+1) = i(k+1) + di;
    end
    
    % update leftmost current from source to node1 (i(1) between source and node1)
    % using: di_source = (dt/L) * (Vsource - v(1) - R * i(1))
    di_s = (dt / L) * (pulse(t) - v(1) - R_eff * i(1) - Rs * i(1));
    i(1) = i(1) + di_s;
    
    % update capacitor voltages: v_k^{n+1} = v_k^{n} + (dt/C)*(i_k - i_{k+1} - G*v_k)
    for k = 1:N
        dv = (dt / C) * (i(k) - i(k+1) - G * v(k));
        v(k) = v(k) + dv;
    end
    
    % apply load at the right end: model as resistor RL between last node and ground
    % current to load approximated via i_load = v(N) / RL; adjust last series current accordingly
    % (we already use currents in updates; optionally damp last current)
    % For numerical stability, enforce last node to leak to ground via RL (a shunt element)
    v(N) = v(N) - (dt / C) * (v(N) / RL);
    
    % record observation nodes
    Vrec(:,n) = v(obs_idx);
end

% ---- Plot results ----
figure('Name','Node voltages','NumberTitle','off');
plot(tvec*1e9, Vrec(1,:), 'LineWidth', 1); hold on;
plot(tvec*1e9, Vrec(2,:), 'LineWidth', 1);
plot(tvec*1e9, Vrec(3,:), 'LineWidth', 1);
plot(tvec*1e9, Vrec(4,:), 'LineWidth', 1);
xlabel('Time (ns)');
ylabel('Voltage (V)');
legend('Node 1 (near source)','Node N/4','Node N/2','Node N (load)');
title('Pulse propagation along lumped LC ladder');
grid on;
