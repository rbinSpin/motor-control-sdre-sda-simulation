%% Updated Main Simulation with SDA-based SDRE Control
clear; clc; close all;

%% 1. Initialize Parameters
param.R = 2.5;
param.Ld = 18e-3;
param.Lq = 32e-3; 
param.lambda_m = 0.30175; 
param.p = 6; 
param.J = 0.00782; 
param.B = 0.02277;

base.Ib = 10; 
base.Vb = 178.82; 
base.wb = 628;
base.Tb = 12.814;

% Calculate PU coefficients (k1 ~ k11)
k = calculate_k_coefficients(param, base);

%% 2. Simulation Settings
dt = 1e-4; % 10kHz sampling
t_end = 0.1; % Reduced to 0.1s for consistency with your HLS test
t_vector = (0:dt:t_end)'; % Ensure column vector
num_steps = length(t_vector);

% Initial states [omega_e_pu; i_q_pu; i_d_pu]
x = [0; 0; 0];
history = zeros(num_steps, 3);
v_history = zeros(num_steps, 2); 
ref_history = zeros(num_steps, 1); 

%% 3. Main Simulation Loop
final_omega_rpm = 900;       
ramp_time = 0.06;            
ramp_slope = final_omega_rpm / ramp_time;

for n = 1:num_steps
    currentTime = t_vector(n);
    
    % --- A. Desired State (Ramp Reference) ---
    if currentTime < ramp_time
        omega_ref_rpm = ramp_slope * currentTime;
    else
        omega_ref_rpm = final_omega_rpm;
    end
    
    ref_history(n) = omega_ref_rpm;
    we_hat = (omega_ref_rpm * pi * param.p / 60) / base.wb; 
    
    % --- B. Controller ---
    u_pu = pmsm_controller(x, we_hat);
    
    % --- C. Plant Dynamics ---
    %TL_pu = 0; 
    %[~, x_next] = ode45(@(t, y) ipmsm_pu_dynamics(t, y, u_pu, k, TL_pu), [0 dt], x);
    %x = x_next(end, :)'; 
    % euler method
    TL_pu = 0; 
    dxdt = ipmsm_pu_dynamics(0, x, u_pu, k, TL_pu);
    x = x + dxdt * dt;
    
    % Save data
    history(n, :) = x';
    v_history(n, :) = u_pu';
    
    if mod(n, 100) == 0
        fprintf('Progress: %.1f%% (%d / %d steps)\n', (n/num_steps)*100, n, num_steps);
    end
end

%% 4. Plotting Results
figure('Color', 'w', 'Name', 'SDRE-SDA Simulation Results');
% Speed
subplot(3,1,1);
act_rpm = history(:,1) * base.wb * 60 / (pi * param.p);
plot(t_vector, ref_history, 'r--', 'LineWidth', 1.2); hold on;
plot(t_vector, act_rpm, 'k-', 'LineWidth', 1.5);
ylabel('Speed (RPM)'); title('SDRE-SDA Control Performance'); 
legend('Ref', 'Act', 'Location', 'southeast'); grid on;

% Currents
subplot(3,1,2);
plot(t_vector, history(:,2), 'r', t_vector, history(:,3), 'b', 'LineWidth', 1.2);
ylabel('Current (pu)'); legend('iq', 'id'); grid on;

% Voltages
subplot(3,1,3);
plot(t_vector, v_history(:,1), 'r', t_vector, v_history(:,2), 'b', 'LineWidth', 1.2);
ylabel('Voltage (pu)'); legend('vq', 'vd'); xlabel('Time (s)'); grid on;

%% 5. Export to CSV (Formatted for HLS Comparison)
% Prepare data for export
Time = t_vector;
Ref_RPM = ref_history;
Act_RPM = act_rpm;
Vq_pu = v_history(:,1);
Vd_pu = v_history(:,2);
Iq_pu = history(:,2); % Convert to Ampere for easier comparison
Id_pu = history(:,3);

% Create table
sim_results_table = table(Time, Ref_RPM, Act_RPM, Vq_pu, Vd_pu, Iq_pu, Id_pu);

% Save to file
output_filename = 'matlab_sim_results.csv';
output_filename2 = '../../data_and_result/sim_result_csv/matlab_sim_results.csv';
writetable(sim_results_table, output_filename);
writetable(sim_results_table, output_filename2);
fprintf('Simulation complete. Data exported to %s\n', output_filename);