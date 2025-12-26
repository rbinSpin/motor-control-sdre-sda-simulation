function u_pu = pmsm_controller(x, we_hat, param, base, k, Q, R)
% This is the top-level function for HDL code generation.
% It encapsulates the entire control algorithm.
%
% Inputs:
%   x       - Current state vector [omega_e_pu; i_q_pu; i_d_pu]
%   we_hat  - Desired electrical speed (pu)
%   param   - Physical parameters of the motor
%   base    - Base values for PU conversion
%   k       - Pre-calculated PU coefficients
%   Q, R    - Weighting matrices for SDRE
%
% Output:
%   u_pu    - Output voltage vector [vq_pu; vd_pu]

% For HDL Coder, all inputs need to be defined. We pass structs like param, 
% base, and k for now, but for hardware, their few used values should be 
% passed as separate constants.

% --- A. Desired State Calculation (based on we_hat) ---
% In a real scenario, Te_target might be a direct input instead of being derived from we_hat.
Te_target = Te_target_calculation(we_hat,k.k2);

[id_hat, iq_hat] = mtpa_calculation(Te_target * base.Tb, param, base);

x_hat = [we_hat; iq_hat; id_hat];

% --- B. Feedforward Voltage (u_hat) ---
[vq_hat, vd_hat] = calculate_feedforward_voltage(k, iq_hat, we_hat, id_hat);
u_hat = [vq_hat; vd_hat];

% --- C. Construct A_bar Matrix (State-Dependent) ---
kk = [k.k1 k.k2 k.k3 k.k4 k.k5 k.k6 k.k7 k.k8 k.k9 k.k10 k.k11];
[A_bar, B] = construct_AB_matrices(kk, x, we_hat, id_hat);

% --- D. SDA Algorithm (Solving ARE) ---
% The sda_riccati_solver is an iterative function. It needs to be HDL-compatible.
P = sda_riccati_solver(A_bar, B, Q, R, 10); % 10 is max iterations

u_pu = calculate_control_output(R,B,P,x,x_hat,u_hat);
end
