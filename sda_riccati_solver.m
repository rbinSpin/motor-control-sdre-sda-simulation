function P = sda_riccati_solver(A_in, B_in, Q_in, R_in, max_iter, T_matrix, FIMATH_SETTINGS)
    % SDA solver for A'P + PA - PBR^-1B'P + Q = 0
    % Fixed-point version
    
    % Attach fimath to all fi objects in this function
    fipref('DataTypeOverride', 'ForceOff');
    globalfimath(FIMATH_SETTINGS);
    
    % --- Cast inputs to fixed-point ---
    A = fi(A_in, T_matrix);
    B = fi(B_in, T_matrix);
    Q = fi(Q_in, T_matrix);
    R = fi(R_in, T_matrix);

    % Step 1: Physical settings and Matrix G
    % G = B * (R \ B'); % This fails for fi objects
    % Manual inversion for diagonal R matrix
    diag_R = diag(R);
    one_fi = fi(1, numerictype(diag_R), fimath(diag_R));
    inv_diag_R = one_fi ./ diag_R;
    R_inv = diag(inv_diag_R);
    G = B * R_inv * B';
    I = fi(eye(size(A)), numerictype(1, 2, 0)); % Identity matrix fi(1)
    I_full = fi(eye(size(A)), T_matrix); % Full precision identity
    
    gamma_val = 500;
    gamma = fi(gamma_val, T_matrix); % Stabilization parameter

    % Step 2: Initialization (Cayley Transform)
    % A_gamma = A - gamma*I 必須能完整表達 -500 左右的數值
    A_gamma = fi(double(A) - gamma * eye(size(A_in)), T_matrix);
    
    % Pre-calculating common terms for efficiency
    % Replace inv() with fixed.qrMatrixSolve for fi support
    X = fixed.qrMatrixSolve(A_gamma', Q);
    W_inner_term = A_gamma - G * X;
    W = fixed.qrMatrixSolve(W_inner_term, I_full);
    
    A0 = I_full + 2 * gamma * W;
    G0 = 2 * gamma * W * G * fixed.qrMatrixSolve(A_gamma', I_full);
    H0 = 2 * gamma * fixed.qrMatrixSolve(A_gamma', Q) * W;
    
    % Step 3: Doubling Iteration Loop
    Ak = A0; Gk = G0; Hk = H0;
    
    for k = 1:max_iter
        % Using the doubling formulas
        inv_I_GH = fixed.qrMatrixSolve(I_full + Gk * Hk, I_full);
        
        Ak_next = Ak * inv_I_GH * Ak;
        Gk_next = Gk + Ak * inv_I_GH * Gk * Ak';
        Hk_next = Hk + Ak' * Hk * inv_I_GH * Ak;
        
        % Update matrices
        Ak = Ak_next;
        Gk = Gk_next;
        Hk = Hk_next;
        
        % Convergence check (using double for norm)
        if norm(double(Ak), 'fro') < 1e-12
            break;
        end
    end
    
    % The solution P is the converged Hk matrix
    P = Hk;
    
    % Reset global fimath
    globalfimath;
end