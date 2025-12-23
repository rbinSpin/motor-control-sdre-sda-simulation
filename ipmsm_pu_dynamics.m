function dxdt = ipmsm_pu_dynamics(~, x_in, u_in, k, TL_pu, T_signal, FIMATH_SETTINGS)
    % Fixed-point version of IPMSM dynamics
    % State vector x = [omega_e_pu; i_q_pu; i_d_pu]
    % Input vector u = [Vq_pu; Vd_pu]
    
    % Attach fimath for this function's scope
    fipref('DataTypeOverride', 'ForceOff');
    globalfimath(FIMATH_SETTINGS);
    
    % --- Cast inputs to fixed-point ---
    x = fi(x_in, T_signal);
    u = fi(u_in, T_signal);
    
    we = x(1); % Electrical angular velocity (PU) 
    iq = x(2); % q-axis current (PU) 
    id = x(3); % d-axis current (PU) 
    
    vq = u(1); % q-axis voltage (PU) 
    vd = u(2); % d-axis voltage (PU) 

    % Mechanical dynamic: d(we)/dt 
    % Te = k.k1*iq + k.k11*id*iq
    % d_we = Te - k.k2*we - k.k3*TL_pu
    d_we = k.k1*iq - k.k2*we - k.k3*TL_pu + k.k11*id*iq;
    
    % q-axis current dynamic: d(iq)/dt 
    d_iq = -k.k4*iq - k.k5*we + k.k6*vq - k.k10*we*id;
    
    % d-axis current dynamic: d(id)/dt 
    d_id = -k.k7*id + k.k8*vd + k.k9*we*iq;

    dxdt = [d_we; d_iq; d_id];
    
    % Reset global fimath
    globalfimath;
end