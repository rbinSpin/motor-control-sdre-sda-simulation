function k = calculate_k_coefficients(param, base)
    % Calculate the PU system coefficients
    
    % Mechanical coefficients
    % k1_bar: Torque constant effect on acceleration 
    k.k1 = (3/2) * (1/param.J) * (param.p^2/4) * param.lambda_m * (1/base.wb);
    
    % k2_bar: Viscous friction effect 
    k.k2 = param.B / param.J;
    
    % k3_bar: Load torque scaling 
    k.k3 = (param.p / (2 * param.J)) * (1 / base.wb);
    
    % Electrical coefficients (q-axis)
    k.k4 = param.R / param.Lq; % Stator resistance effect 
    k.k5 = (param.lambda_m / param.Lq) * (base.wb / base.Ib); % Back-EMF effect 
    k.k6 = (1 / param.Lq) * (base.Vb / base.Ib); % Voltage input gain 
    
    % Electrical coefficients (d-axis)
    k.k7 = param.R / param.Ld; % Stator resistance effect 
    k.k8 = (1 / param.Ld) * (base.Vb / base.Ib); % Voltage input gain 
    
    % Coupling coefficients
    k.k9 = (param.Lq / param.Ld) * base.wb; % Cross-coupling Lq/Ld 
    k.k10 = (param.Ld / param.Lq) * base.wb; % Cross-coupling Ld/Lq 
    
    % Salience torque coefficient (for IPMSM) 
    k.k11 = (3/2) * (1/param.J) * (param.p^2/4) * (param.Ld - param.Lq) * (base.Ib^2 / base.wb);
end