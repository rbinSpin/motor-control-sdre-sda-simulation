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
t_end = 0.1;
t_vector = 0:dt:t_end;
num_steps = length(t_vector);

% Initial states [omega_e_pu; i_q_pu; i_d_pu]
x = [0; 0; 0];
history = zeros(num_steps, 3);
v_history = zeros(num_steps, 2); % To monitor vd, vq

% SDRE Weights
Q = diag([50, 10, 10]); % [we_err, iq_err, id_err]
R = diag([30, 30]);  % [vq, vd]

%% 3. Main Simulation Loop
% 設定斜坡參數
final_omega_rpm = 900;       % 最終目標轉速
ramp_time = 0.06;            % 預計在 ramp_time 秒內升到目標轉速
ramp_slope = final_omega_rpm / ramp_time;

for n = 1:num_steps
    currentTime = t_vector(n);
    
    % --- A. Desired State (Ramp Reference) ---
    % 計算當前的需求轉速
    if currentTime < ramp_time
        omega_ref_rpm = ramp_slope * currentTime;
    else
        omega_ref_rpm = final_omega_rpm;
    end
    
    % 將 RPM 轉換為標么值
    we_hat = (omega_ref_rpm * pi * param.p / 60) / base.wb; 
    
    % 利用物理穩態平衡反推需要的轉矩
    % Note: The reference generation is part of the testbench.
    % In a real system, we might get we_hat from a higher-level controller.
    
    % --- Call Controller (The DUT) ---
    u_pu = pmsm_controller(x, we_hat, param, base, k, Q, R);

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