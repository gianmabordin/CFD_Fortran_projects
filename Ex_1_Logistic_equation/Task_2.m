% =============================================
% PLOT TASK 2: ERROR ANALYSIS AND FITTING IN DP
% =============================================
clear; clc; close all;

% 1. Data loading
data_err = load('task_2_errors.txt');

% 2. Column selection
dt        = data_err(:, 1); 
err_euler = data_err(:, 2);
err_cn    = data_err(:, 3);

% 2. Mathematical fitting
% Linear regression (least squares) performed with the "\" operator
A1 = dt \ err_euler;          
A2 = (dt.^2) \ err_cn;        

% Fitting curves creation
fit_euler = A1 * (dt.^1);
fit_cn    = A2 * (dt.^2);

% 3. Log-Log plot
figure('Name', 'Task 2 - Error Analysis dp', 'Color', 'w');

% Fitting lines
loglog(dt, fit_euler, 'b--', 'LineWidth', 1); hold on;
loglog(dt, fit_cn,    'r--',  'LineWidth', 1);

% Simulation data points
loglog(dt, err_euler, 'b*', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b');
loglog(dt, err_cn,    'r*', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'r');

% 4. Labels
xlabel('Time Step \Delta t [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Absolute Error ||e_N||', 'FontSize', 12, 'FontWeight', 'bold');


% 5. Legend (with coefficient inclusion)
leg_euler_fit = sprintf('Implicit Euler Fit: %.3f \\cdot \\Delta t', A1);
leg_cn_fit    = sprintf('Crank-Nicolson Fit: %.3f \\cdot \\Delta t^2', A2);

legend(leg_euler_fit, leg_cn_fit, 'Implicit Euler (Data)', 'Crank-Nicolson (Data)', ...
       'Location', 'best', 'FontSize', 11);


grid on;
grid minor;
box on;