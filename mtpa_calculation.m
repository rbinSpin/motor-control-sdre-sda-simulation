function [id_ref, iq_ref] = mtpa_calculation(Te_ref, param, base)
    % Calculate MTPA reference currents using a Look-Up Table (LUT).
    % This is the hardware-friendly approach to replace the 'roots' function.
    % The 'param' argument is unused but kept for interface compatibility.

    % A pre-calculated LUT for MTPA, based on the original polynomial. 
    % In a real project, this would be generated offline with higher resolution.
    % Columns: [Torque (Nm), id_physical (A), iq_physical (A)]
    persistent mtpa_lut;
    if isempty(mtpa_lut)
        mtpa_lut = [
             0.0,  0.0,      0.0;
             1.0, -0.025,   0.736;
             2.0, -0.10,    1.47;
             3.0, -0.22,    2.20;
             4.0, -0.40,    2.92;
             5.0, -0.62,    3.65;
             6.0, -0.90,    4.36;
             7.0, -1.22,    5.08;
             8.0, -1.60,    5.78;
             9.0, -2.02,    6.48;
            10.0, -2.50,    7.18
        ];
    end
    
    % For HDL synthesis, it's better to separate columns into individual arrays.
    Te_col = mtpa_lut(:, 1);
    id_col = mtpa_lut(:, 2);
    iq_col = mtpa_lut(:, 3);
    
    % Find the interval in the LUT for the given Te_ref.
    % This loop structure is synthesizable.
    % It finds the index of the lower bound for the interpolation.
    index_low = 1;
    for i = 1:size(mtpa_lut, 1)
        if Te_ref >= Te_col(i)
            index_low = i;
        end
    end
    
    % To prevent out-of-bounds access, if index_low is the last entry,
    % we cannot interpolate, so we set the upper index to be the same.
    if index_low == size(mtpa_lut, 1)
        index_high = index_low;
    else
        index_high = index_low + 1;
    end

    % Get the four points for linear interpolation
    Te_low  = Te_col(index_low);
    id_low  = id_col(index_low);
    iq_low  = iq_col(index_low);
    
    Te_high = Te_col(index_high);
    id_high = id_col(index_high);
    iq_high = iq_col(index_high);

    % Perform linear interpolation
    % Avoid division by zero if Te_high is very close to Te_low
    delta_Te = Te_high - Te_low;
    if delta_Te < 1e-9
        fraction = 0;
    else
        fraction = (Te_ref - Te_low) / delta_Te;
    end
    
    id_ref_physical = id_low + fraction * (id_high - id_low);
    iq_ref_physical = iq_low + fraction * (iq_high - iq_low);

    % Convert physical values to Per-unit (PU) system
    iq_ref = iq_ref_physical / base.Ib;
    id_ref = id_ref_physical / base.Ib;
end