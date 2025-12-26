function P = sda_riccati_solver(A, B, Q, R, max_iter)
    % SDA solver for A'P + PA - PBR^-1B'P + Q = 0
    % Based on SDA_FPGA計算步驟.pptx formulas
    
    % Step 1: Physical settings and Matrix G
    % Replace non-synthesizable matrix division R \ B'
    inv_R = diag([1/R(1,1), 1/R(2,2)]);
    G = B * inv_R * B';
    
    I = eye(size(A));
    gamma = 500; % Stabilization parameter (Cayley transform)
    
    % Step 2: Initialization (Cayley Transform)
    % As seen in the PPTX: Initializing A0, G0, H0
    A_gamma = A - gamma * I;
    
    % Pre-calculating common terms for efficiency using HDL-compatible QR decomposition
    inv_A_gamma_t = inv3x3_qr_hdl(A_gamma');
    
    W_interim = A_gamma - G * (inv_A_gamma_t * Q);
    W = inv3x3_qr_hdl(W_interim);
    
    A0 = I + 2 * gamma * W;
    G0 = 2 * gamma * W * G * inv_A_gamma_t;
    H0 = 2 * gamma * (inv_A_gamma_t * Q) * W;
    
    % Step 3: Doubling Iteration Loop
    Ak = A0; Gk = G0; Hk = H0;
    
    % Use a persistent flag for convergence to be HDL-friendly
    converged = false;
    
    for k = 1:max_iter
        % Only perform calculations if not yet converged
        if ~converged
            % Using the doubling formulas from your PPTX, with HDL-compatible inversion
            inv_I_GH = inv3x3_qr_hdl(I + Gk * Hk);
            
            Ak_next = Ak * inv_I_GH * Ak;
            Gk_next = Gk + Ak * inv_I_GH * Gk * Ak';
            Hk_next = Hk + Ak' * Hk * inv_I_GH * Ak;
            
            % Update matrices
            Ak = Ak_next;
            Gk = Gk_next;
            Hk = Hk_next;
            
            % Convergence check (Hardware-friendly version)
            % Sum of absolute values instead of Frobenius norm
            if sum(abs(Ak(:))) < 1e-10
                converged = true;
            end
        end
    end
    
    % The solution P is the converged Hk matrix
    P = Hk;
end