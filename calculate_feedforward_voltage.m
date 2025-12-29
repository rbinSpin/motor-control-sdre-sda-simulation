function [vq_hat, vd_hat] = calculate_feedforward_voltage(k, iq_hat, we_hat, id_hat)
% CALCULATE_FEEDFORWARD_VOLTAGE Calculates the feedforward voltage components.
%   [vq_hat, vd_hat] = calculate_feedforward_voltage(k, iq_hat, we_hat, id_hat)
%   Calculates the q-axis (vq_hat) and d-axis (vd_hat) feedforward voltage
%   components based on motor coefficients (k) and desired current/speed.

vq_hat = (k(4)*iq_hat + k(5)*we_hat + k(10)*we_hat*id_hat)/k(6);
vd_hat = (k(7)*id_hat - k(9)*we_hat*iq_hat)/k(8);

end
