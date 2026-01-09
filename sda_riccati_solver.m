function P = sda_riccati_solver(A, B, Q, R, max_iter, gamma)
    % SDA solver for A'P + PA - PBR^-1B'P + Q = 0
    % Based on SDA_FPGA計算步驟.pptx formulas
    
    % Step 1: Physical settings and Matrix G
    % Replace non-synthesizable matrix division R \ B'
    inv_R = diag([1/R(1,1), 1/R(2,2)]);
    G = B * inv_R * B';
    
    I = eye(size(A));
    
    % Step 2: Initialization (Cayley Transform)
    % As seen in the PPTX: Initializing A0, G0, H0
    A_gamma = A - gamma * I;
    
    % Pre-calculating common terms for efficiency using HDL-compatible LU decomposition
    inv_A_gamma_t = inv3x3_lu_alaff_hdl(A_gamma');
    blue = A_gamma' + Q*inv_A_gamma_t*G;
    green = Q*inv_A_gamma_t;
    red = G*inv3x3_lu_alaff_hdl(blue);
    inv_A_gamma = inv3x3_lu_alaff_hdl(A_gamma);
    
    A0 = I + 2 * gamma * inv3x3_lu_alaff_hdl(blue);
    G0 = 2 * gamma * inv_A_gamma * red;
    H0 = 2 * gamma * inv3x3_lu_alaff_hdl(blue) * green;
    
    % Step 3: Doubling Iteration Loop
    Ak = A0; Gk = G0; Hk = H0;
    
    for k = 1:max_iter
        % Using the doubling formulas from your PPTX, with HDL-compatible inversion
        inv_I_GH = inv3x3_lu_alaff_hdl(I + Gk * Hk);
        inv_I_HG = inv3x3_lu_alaff_hdl(I + Hk * Gk);
        
        Ak_next = Ak * inv_I_GH * Ak;
        Gk_next = Gk + Ak * Gk * inv_I_HG * Ak';
        Hk_next = Hk + Ak' * inv_I_HG * Hk * Ak;
        
        % Update matrices
        Ak = Ak_next;
        Gk = Gk_next;
        Hk = Hk_next;
    end
    
    % The solution P is the converged Hk matrix
    P = Hk;
end