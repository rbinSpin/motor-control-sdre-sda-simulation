function [id_ref, iq_ref] = mtpa_calculation(Te_ref, param, base)
    % Calculate MTPA reference currents based on meeting0627.pdf (p.11)
    % and meeting0620.pdf (p.23)
    
    % Get physical parameters for calculation
    Ld = param.Ld;
    Lq = param.Lq;
    lambda_m = param.lambda_m;
    P = param.p / 2; % Pole pairs
    Ldq = Ld - Lq;

    % Setup the cubic equation coefficients for iq:
    % Equation: 1.5 * P * (Ld - Lq)^2 / lambda_m * iq^3 + 1.5 * P * lambda_m * iq - Te = 0
    % [cite: 588, 592, 672-686]
    
    coeff_p1 = 1.5 * P * (Ldq^2) / lambda_m;
    coeff_p2 = 0;
    coeff_p3 = 1.5 * P * lambda_m;
    coeff_p4 = -Te_ref;
    
    % Solve for real roots of the cubic equation
    poly_coeffs = [coeff_p1, coeff_p2, coeff_p3, coeff_p4];
    r = roots(poly_coeffs);
    
    % Filter out imaginary solutions and pick the maximum real one [cite: 696-697]
    real_solutions = r(imag(r) == 0);
    iq_ref_physical = max(real_solutions);
    
    % Calculate id_ref using the derived MTPA relationship:
    % id = (Ld - Lq) / lambda_m * iq^2 [cite: 588, 592, 700]
    id_ref_physical = (Ldq / lambda_m) * (iq_ref_physical^2);
    
    % Convert physical values to Per-unit (PU) system
    iq_ref = iq_ref_physical / base.Ib;
    id_ref = id_ref_physical / base.Ib;
end