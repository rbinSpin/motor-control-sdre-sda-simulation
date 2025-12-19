%% Updated Main Simulation with SDA-based SDRE Control
clear; clc; close all;

%% 1. Initialize Parameters
param.R = 2.5; param.Ld = 18e-3; param.Lq = 32e-3; 
param.lambda_m = 0.30175; param.p = 6; param.J = 0.00782; param.B = 0.02277;

base.Ib = 10; base.Vb = 178.82; base.wb = 628;
base.Tb = 12.814;

% Calculate PU coefficients (k1 ~ k11)
k = calculate_k_coefficients(param, base);

%% 2. Simulation Settings
dt = 1e-4; % 10kHz sampling
t_end = 0.5;
t_vector = 0:dt:t_end;
num_steps = length(t_vector);

% Initial states [omega_e_pu; i_q_pu; i_d_pu]
x = [0; 0; 0]; 
history = zeros(num_steps, 3);
v_history = zeros(num_steps, 2); % To monitor vd, vq

% SDRE Weights
Q = diag([50, 10, 10]); % [we_err, iq_err, id_err]
R = diag([1, 1]);           % [vq, vd]

%% 3. Main Simulation Loop
for n = 1:num_steps
    % --- A. Desired State (Constant Reference) ---
    Te_target = 2.0; 
    omega_ref_rpm = 900;
    
    % Get desired currents from MTPA
    [id_hat, iq_hat] = mtpa_calculation(Te_target, param, base);
    we_hat = (omega_ref_rpm * pi * param.p / 60) / base.wb; % Convert RPM to PU
    x_hat = [we_hat; iq_hat; id_hat];
    
    % --- B. Feedforward Voltage (u_hat) ---
    vq_hat = k.k4*iq_hat + k.k5*we_hat + k.k10*we_hat*id_hat;
    vd_hat = k.k7*id_hat - k.k9*we_hat*iq_hat;
    u_hat = [vq_hat; vd_hat];
    
    % --- C. Construct A_bar Matrix (3x3) ---
    A_bar = [ -k.k2,                 (k.k1 + k.k11*id_hat),   k.k11*x(2);
              -(k.k5 + k.k10*x(3)),  -k.k4,                   -k.k10*we_hat;
              k.k9*x(2),             k.k9*we_hat,             -k.k7 ];
              
    B = [ 0,    0;
          k.k6, 0;
          0,    k.k8 ];
    
    % --- D. SDA Algorithm (Solving ARE) ---
    P = sda_riccati_solver(A_bar, B, Q, R, 10);
    K = (R \ B') * P;
    
    % --- E. Calculate Feedback & Total Input ---
    x_tilde = x - x_hat;
    u_tilde = -K * x_tilde;

    % 計算總輸出電壓
    u_pu = u_hat + u_tilde;
    
    % 加上 Saturation 限制 (標么值限制在 [-1.1, 1.1])
    limit = 1.1;
    u_pu(1) = max(min(u_pu(1), limit), -limit); % vq 限制
    u_pu(2) = max(min(u_pu(2), limit), -limit); % vd 限制
    
    % --- F. Plant Dynamics (Motor Model) ---
    TL_pu = 0; % Assuming constant load torque
    [~, x_next] = ode45(@(t, y) ipmsm_pu_dynamics(t, y, u_pu, k, TL_pu), [0 dt], x);
    x = x_next(end, :)'; 
    
    % 歐拉法
    %dxdt = ipmsm_pu_dynamics(0, x, u_pu, k, TL_pu);
    %x = x + dxdt * dt;
    
    % Save data
    history(n, :) = x';
    v_history(n, :) = u_pu';

    % 每 500 步顯示一次進度
    if mod(n, 5) == 0
        fprintf('目前進度: %.1f%% (第 %d 步 / 共 %d 步)\n', (n/num_steps)*100, n, num_steps);
    end
end

%% 4. Plotting Results
figure('Color', 'w');
subplot(3,1,1);
plot(t_vector, history(:,1) * base.wb * 60 / (pi * param.p), 'k', 'LineWidth', 1.5);
ylabel('Speed (RPM)'); title('SDRE-SDA Control Performance'); grid on;

subplot(3,1,2);
plot(t_vector, history(:,2), 'r', t_vector, history(:,3), 'b');
ylabel('Current (pu)'); legend('iq', 'id'); grid on;

subplot(3,1,3);
plot(t_vector, v_history(:,1), 'r', t_vector, v_history(:,2), 'b');
ylabel('Voltage (pu)'); legend('vq', 'vd'); xlabel('Time (s)'); grid on;