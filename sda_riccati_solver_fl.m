function P = sda_riccati_solver_fl(A, B, Q, R, max_iter)
    % SDA solver for A'P + PA - PBR^-1B'P + Q = 0
    % Based on SDA_FPGA計算步驟.pptx formulas
    
    % Step 1: Physical settings and Matrix G
    G = B * (R \ B');
    I = eye(size(A));
    gamma = 500; % Stabilization parameter (Cayley transform)
    
    % Step 2: Initialization (Cayley Transform)
    % As seen in the PPTX: Initializing A0, G0, H0
    A_gamma = A - gamma * I;
    
    % Pre-calculating common terms for efficiency
    % W = (A_gamma - G * (A_gamma' \ Q)) \ I
    W = (A_gamma - G * (A_gamma' \ Q)) \ I;
    
    A0 = I + 2 * gamma * W;
    G0 = 2 * gamma * W * G * (A_gamma' \ I);
    H0 = 2 * gamma * (A_gamma' \ Q) * W;
    
    % Step 3: Doubling Iteration Loop
    Ak = A0; Gk = G0; Hk = H0;
    
    for k = 1:max_iter
        % Using the doubling formulas from your PPTX
        inv_I_GH = (I + Gk * Hk) \ I;
        
        Ak_next = Ak * inv_I_GH * Ak;
        Gk_next = Gk + Ak * inv_I_GH * Gk * Ak';
        Hk_next = Hk + Ak' * Hk * inv_I_GH * Ak;
        
        % Update matrices
        Ak = Ak_next;
        Gk = Gk_next;
        Hk = Hk_next;
        
        % Convergence check
        if norm(Ak, 'fro') < 1e-12
            break;
        end
    end
    
    % The solution P is the converged Hk matrix
    P = Hk;
end