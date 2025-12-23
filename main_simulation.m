%% Updated Main Simulation with SDA-based SDRE Control
clear; clc; close all;

%% 1. Initialize Parameters
param.R = 2.5; param.Ld = 18e-3; param.Lq = 32e-3; 
param.lambda_m = 0.30175; param.p = 6; param.J = 0.00782; param.B = 0.02277;

base.Ib = 10; base.Vb = 178.82; base.wb = 628;
base.Tb = 12.814;

% Calculate PU coefficients (k1 ~ k11)
k = calculate_k_coefficients(param, base);

%% 2. Fixed-Point Settings
% Define data types for signals and calculations
% T_signal: For per-unit signals. Increased from 16 to 32 bits for higher precision.
T_signal = numerictype('Signed', true, 'WordLength', 32, 'FractionLength', 28);

% T_matrix: For matrices. Increased from 32 to 48 bits.
T_matrix = numerictype('Signed', true, 'WordLength', 48, 'FractionLength', 40);

% T_prod: For products and accumulations. Increased from 32 to 48 bits.
T_prod = numerictype('Signed', true, 'WordLength', 48, 'FractionLength', 44);

% Associate these types with a fimath object to control rounding and overflow behavior.
% Rounding to nearest, saturate on overflow.
FIMATH_SETTINGS = fimath('RoundingMethod', 'Floor', 'OverflowAction', 'Saturate', ...
    'ProductMode', 'SpecifyPrecision', 'ProductWordLength', T_prod.WordLength, 'ProductFractionLength', T_prod.FractionLength, ...
    'SumMode', 'SpecifyPrecision', 'SumWordLength', T_prod.WordLength, 'SumFractionLength', T_prod.FractionLength);

%% Load MTPA Look-up Table
load('mtpa_lut.mat');

%% 2. Simulation Settings
dt = 1e-4; % 10kHz sampling
t_end = 0.005;
t_vector = 0:dt:t_end;
num_steps = length(t_vector);

% --- Floating-Point (FL) Initialization ---
x_fl = [0; 0; 0];
history_fl = zeros(num_steps, 3);
v_history_fl = zeros(num_steps, 2);

% --- Fixed-Point (FX) Initialization ---
fipref('DataTypeOverride', 'ForceOff');
globalfimath(FIMATH_SETTINGS);
x_fx = fi([0; 0; 0], T_signal); 
history_fx = zeros(num_steps, 3);
v_history_fx = zeros(num_steps, 2); 

% SDRE Weights
Q_fl = diag([5000, 10, 10]); 
R_fl = diag([30, 30]);

%% 3. Main Simulation Loop
fprintf('Starting simulation...\n');
% 設定斜坡參數
final_omega_rpm = 900;
ramp_time = 0.006;
ramp_slope = final_omega_rpm / ramp_time;

for n = 1:num_steps
    currentTime = t_vector(n);
    
    % --- A. Desired State (Common for both FL and FX) ---
    if currentTime < ramp_time
        omega_ref_rpm = ramp_slope * currentTime;
    else
        omega_ref_rpm = final_omega_rpm;
    end
    
    we_hat_fl = (omega_ref_rpm * pi * param.p / 60) / base.wb; 
    Te_target_fl = k.k2 * we_hat_fl;
    
    iq_hat_fl = interp1(Te_lut, iq_lut, Te_target_fl * base.Tb, 'linear', 'extrap');
    id_hat_fl = interp1(Te_lut, id_lut, Te_target_fl * base.Tb, 'linear', 'extrap');
    
    x_hat_fl = [we_hat_fl; iq_hat_fl; id_hat_fl];
    
    % --- B. Feedforward Voltage (Common for both FL and FX) ---
    vq_hat_fl = (k.k4*iq_hat_fl + k.k5*we_hat_fl + k.k10*we_hat_fl*id_hat_fl)/k.k6;
    vd_hat_fl = (k.k7*id_hat_fl - k.k9*we_hat_fl*iq_hat_fl)/k.k8;
    u_hat_fl = [vq_hat_fl; vd_hat_fl];

    % --- C. Construct A_bar Matrix (State-dependent) ---
    A_bar_fl = [ -k.k2, (k.k1 + k.k11*id_hat_fl), k.k11*x_fl(2);
                 -(k.k5 + k.k10*x_fl(3)), -k.k4, -k.k10*we_hat_fl;
                  k.k9*x_fl(2), k.k9*we_hat_fl, -k.k7 ];
    B_fl = [ 0, 0; k.k6, 0; 0, k.k8 ];
    
    % =====================================================================
    % --- FLOATING-POINT (FL) PATH ---
    % =====================================================================
    P_fl = sda_riccati_solver_fl(A_bar_fl, B_fl, Q_fl, R_fl, 10);
    K_fl = (R_fl \ B_fl') * P_fl;
    x_tilde_fl = x_fl - x_hat_fl;
    u_tilde_fl = -K_fl * x_tilde_fl;
    u_pu_fl = u_hat_fl + u_tilde_fl;
    u_pu_fl(1) = max(min(u_pu_fl(1), 1.1), -1.1);
    u_pu_fl(2) = max(min(u_pu_fl(2), 1.1), -1.1);
    dxdt_fl = ipmsm_pu_dynamics_fl(0, x_fl, u_pu_fl, k, 0);
    x_fl = x_fl + dxdt_fl * dt;
    history_fl(n, :) = x_fl';
    v_history_fl(n, :) = u_pu_fl';

    % =====================================================================
    % --- FIXED-POINT (FX) PATH ---
    % =====================================================================
    x_hat_fx = fi(x_hat_fl, T_signal);
    u_hat_fx = fi(u_hat_fl, T_signal);
    
    A_bar_fx_fl = [ -k.k2, (k.k1 + k.k11*id_hat_fl), k.k11*double(x_fx(2));
                    -(k.k5 + k.k10*double(x_fx(3))), -k.k4, -k.k10*we_hat_fl;
                     k.k9*double(x_fx(2)), k.k9*we_hat_fl, -k.k7 ];

    A_bar_fx = fi(A_bar_fx_fl, T_matrix);
    B_fx = fi(B_fl, T_matrix);
    Q_fx = fi(Q_fl, T_matrix);
    R_fx = fi(R_fl, T_matrix);
    
    P_fx = sda_riccati_solver(A_bar_fx, B_fx, Q_fx, R_fx, 10, T_matrix, FIMATH_SETTINGS);
    
    % K calculation must also be fixed-point, with manual inversion
    diag_R_fx = diag(R_fx);
    one_fi_K = fi(1, numerictype(diag_R_fx), fimath(diag_R_fx));
    inv_diag_R_fx = one_fi_K ./ diag_R_fx;
    R_inv_fx = diag(inv_diag_R_fx);
    K_fx = R_inv_fx * B_fx' * P_fx;
    x_tilde_fx = x_fx - x_hat_fx;
    u_tilde_fx = -K_fx * x_tilde_fx;
    u_pu_fx = u_hat_fx + u_tilde_fx;
    
    limit_fx = fi(1.1, T_signal);
    u_pu_fx(1) = max(min(u_pu_fx(1), limit_fx), -limit_fx);
    u_pu_fx(2) = max(min(u_pu_fx(2), limit_fx), -limit_fx);
    
    dxdt_fx = ipmsm_pu_dynamics(0, x_fx, u_pu_fx, k, 0, T_signal, FIMATH_SETTINGS);
    x_fx = x_fx + dxdt_fx * dt;
    
    history_fx(n, :) = double(x_fx');
    v_history_fx(n, :) = double(u_pu_fx');

    if mod(n, 10) == 0
        fprintf('Progress: %.1f%% (Step %d / %d)\n', (n/num_steps)*100, n, num_steps);
    end
end
fprintf('Simulation finished.\n');
% Clear the global fimath settings
globalfimath;

%% Save Simulation Data for Analysis
save('simulation_results.mat', 't_vector', 'history_fl', 'history_fx', 'v_history_fl', 'v_history_fx');
fprintf('Simulation results saved to simulation_results.mat for analysis.\n');

%% 4. Plotting Results
figure('Color', 'w', 'Name', 'Fixed-Point vs Floating-Point Comparison');

% ----- Speed Plot -----
ax1 = subplot(3,3,1:3);
plot(t_vector, history_fl(:,1) * base.wb * 60 / (pi * param.p), 'b:', 'LineWidth', 2.5); hold on;
plot(t_vector, history_fx(:,1) * base.wb * 60 / (pi * param.p), 'k');
grid on; title('Motor Speed'); ylabel('Speed (RPM)');
legend('Float', 'Fixed', 'Location', 'southeast');

% ----- q-axis Current Plot -----
ax2 = subplot(3,3,4);
plot(t_vector, history_fl(:,2), 'b:', 'LineWidth', 2.5); hold on;
plot(t_vector, history_fx(:,2), 'r');
grid on; title('q-axis Current'); ylabel('iq (pu)');

% ----- d-axis Current Plot -----
ax3 = subplot(3,3,5);
plot(t_vector, history_fl(:,3), 'b:', 'LineWidth', 2.5); hold on;
plot(t_vector, history_fx(:,3), 'r');
grid on; title('d-axis Current'); ylabel('id (pu)');

% ----- Error Plot -----
ax4 = subplot(3,3,6);
error_rpm = (history_fl(:,1) - history_fx(:,1)) * base.wb * 60 / (pi * param.p);
plot(t_vector, error_rpm, 'm');
grid on; title('Speed Error'); ylabel('Error (RPM)');

% ----- Voltage Plots -----
ax5 = subplot(3,3,7:8);
plot(t_vector, v_history_fl(:,1), 'b:', 'LineWidth', 2.5); hold on;
plot(t_vector, v_history_fx(:,1), 'r');
plot(t_vector, v_history_fl(:,2), 'c:', 'LineWidth', 2.5);
plot(t_vector, v_history_fx(:,2), 'k');
grid on; title('Voltages'); ylabel('Voltage (pu)'); xlabel('Time (s)');
legend('Vq Float', 'Vq Fixed', 'Vd Float', 'Vd Fixed');

linkaxes([ax1, ax2, ax3, ax4, ax5], 'x');
sgtitle('Floating-Point vs. Fixed-Point Simulation Results');