function u_total = sda_control_logic(x, x_hat, k, Q, R)
    % x: [we; iq; id] current PU
    % x_hat: [we_hat; iq_hat; id_hat] desired PU (Constant)
    
    % 1. Calculate Feedforward Voltage u_hat (Desired state maintenance)
    vq_hat = k.k4*x_hat(2) + k.k5*x_hat(1) + k.k10*x_hat(1)*x_hat(3);
    vd_hat = k.k7*x_hat(3) - k.k9*x_hat(1)*x_hat(2);
    u_hat = [vq_hat; vd_hat];
    
    % 2. Construct your refined 3x3 bar(A) matrix
    A_bar = [ -k.k2,            (k.k1 + k.k11*x_hat(3)),  k.k11*x(2);
              -(k.k5 + k.k10*x(3)), -k.k4,                -k.k10*x_hat(1);
              (k.k9*x(2) + k.k9*x_hat(1)), k.k9*x_hat(1), -k.k7 ];
              
    % 3. Define B matrix (3x2)
    B = [ 0,    0;
          k.k6, 0;
          0,    k.k8 ];
    
    % 4. Solve ARE using SDA algorithm
    % Typically 5-10 iterations are enough for convergence
    P = sda_riccati_solver(A_bar, B, Q, R, 10);
    
    % 5. Calculate Feedback Gain K and Error Compensation u_tilde
    K = (R \ B') * P;
    x_tilde = x - x_hat;
    u_tilde = -K * x_tilde;
    
    % 6. Final Control Output
    u_total = u_hat + u_tilde;
end