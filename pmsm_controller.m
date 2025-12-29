function u_pu = pmsm_controller(x, we_hat)
% This is the top-level function for HDL code generation.
% It encapsulates the entire control algorithm.
%
% Inputs: all inputs are arrays.
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

persistent base_array k_array Q R
    
    if isempty(base_array)
        base_array = [10.0000  178.8200  628.0000   12.8140];
    end

    if isempty(k_array)
        k_array = [0.8295  2.9118  0.6109  78.1250  592.1844  558.8125  138.8889  993.4444  1116.4444  353.2500  -3.8485];
    end

    if isempty(Q)
        Q = diag([5000, 10, 10]);
    end

    if isempty(R)
        R = diag([30, 30]);
    end
% --- A. Desired State Calculation (based on we_hat) ---
% In a real scenario, Te_target might be a direct input instead of being derived from we_hat.
Te_target = Te_target_calculation(we_hat,k_array(2));

[id_hat, iq_hat] = mtpa_calculation(Te_target * base_array(4), base_array);

x_hat = [we_hat; iq_hat; id_hat];

% --- B. Feedforward Voltage (u_hat) ---
[vq_hat, vd_hat] = calculate_feedforward_voltage(k_array, iq_hat, we_hat, id_hat);
u_hat = [vq_hat; vd_hat];

% --- C. Construct A_bar Matrix (State-Dependent) ---
[A_bar, B] = construct_AB_matrices(k_array, x, we_hat, id_hat);

% --- D. SDA Algorithm (Solving ARE) ---
% The sda_riccati_solver is an iterative function. It needs to be HDL-compatible.
P = sda_riccati_solver(A_bar, B, Q, R, 10); % 10 is max iterations

u_pu = calculate_control_output(R,B,P,x,x_hat,u_hat);
end
