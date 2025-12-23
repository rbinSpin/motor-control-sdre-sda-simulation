%% Generate MTPA Look-Up Table

clear; clc; close all;

% This script generates a look-up table (LUT) for the MTPA (Maximum Torque
% Per Ampere) calculation. It pre-computes the optimal id and iq currents
% for a range of torque demands and saves the results into a .mat file.

%% 1. Initialize Motor Parameters (must match main_simulation)
param.R = 2.5; param.Ld = 18e-3; param.Lq = 32e-3; 
param.lambda_m = 0.30175; param.p = 6; param.J = 0.00782; param.B = 0.02277;

base.Ib = 10; base.Vb = 178.82; base.wb = 628;
base.Tb = 12.814;

%% 2. Generate LUT Data
% Define the torque vector for the LUT.
% Let's span a range slightly larger than the base torque.
num_points = 256;
Te_lut = linspace(0, 1.2 * base.Tb, num_points); % Torque in Nm

% Pre-allocate arrays for the results
id_lut = zeros(1, num_points);
iq_lut = zeros(1, num_points);

% Loop through each torque value and calculate the optimal currents
for i = 1:num_points
    [id_lut(i), iq_lut(i)] = mtpa_calculation(Te_lut(i), param, base);
end

%% 3. Save LUT to .mat file
save('mtpa_lut.mat', 'Te_lut', 'id_lut', 'iq_lut');

fprintf('MTPA look-up table generated and saved to mtpa_lut.mat\n');

%% 4. Plot the LUT for verification
figure('Color', 'w', 'Name', 'MTPA Look-Up Table');
subplot(2,1,1);
plot(Te_lut, id_lut, 'b', 'LineWidth', 1.5);
title('MTPA Look-Up Table');
ylabel('id_ref (A)');
grid on;

subplot(2,1,2);
plot(Te_lut, iq_lut, 'r', 'LineWidth', 1.5);
ylabel('iq_ref (A)');
xlabel('Torque Command (Nm)');
grid on;

sgtitle('Generated MTPA id/iq Reference Currents');
