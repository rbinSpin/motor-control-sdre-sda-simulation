function M_inv = inv3x3_lu_alaff_hdl(M)
% inv3x3_lu_alaff_hdl Hardware-friendly 3x3 matrix inversion using LU decomposition.
%   This function solves M*x = I for x, where x is the inverse of M.
%   The procedure is:
%   1. Decompose M into L*U (Doolittle method).
%   2. Solve L*y = I for y using forward substitution.
%   3. Solve U*x = y for x using backward substitution.

    L = eye(3);
    U = zeros(3);
    
    % --- Step 1: LU Decomposition (Doolittle Algorithm) ---
    % This part can be sensitive to zeros on the diagonal.
    % Adding eps for numerical stability.
    
    % First row of U
    U(1,1) = M(1,1) + eps;
    U(1,2) = M(1,2);
    U(1,3) = M(1,3);

    % First column of L
    L(2,1) = M(2,1) / U(1,1);
    L(3,1) = M(3,1) / U(1,1);
    
    % Second row of U
    U(2,2) = M(2,2) - L(2,1) * U(1,2) + eps;
    U(2,3) = M(2,3) - L(2,1) * U(1,3);

    % Second column of L
    L(3,2) = (M(3,2) - L(3,1) * U(1,2)) / U(2,2);
    
    % Third row of U
    U(3,3) = M(3,3) - L(3,1) * U(1,3) - L(3,2) * U(2,3) + eps;

    % --- Steps 2 & 3: Forward/Backward Substitution for each column of I ---
    I = eye(3);
    M_inv = zeros(3);
    
    for i = 1:3
        b = I(:, i);
        
        % Step 2: Solve Ly = b for y (Forward Substitution)
        y = zeros(3,1);
        y(1) = b(1);
        y(2) = b(2) - L(2,1) * y(1);
        y(3) = b(3) - L(3,1) * y(1) - L(3,2) * y(2);
        
        % Step 3: Solve Ux = y for x (Backward Substitution)
        x = zeros(3,1);
        x(3) = y(3) / U(3,3);
        x(2) = (y(2) - U(2,3) * x(3)) / U(2,2);
        x(1) = (y(1) - U(1,2) * x(2) - U(1,3) * x(3)) / U(1,1);
        
        M_inv(:, i) = x;
    end
end
