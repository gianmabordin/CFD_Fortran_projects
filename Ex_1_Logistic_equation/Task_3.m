% ============================================
% PLOT TASK 3: ERROR ANALYSIS SINGLE PRECISION
% ============================================
clear; clc; close all;

% 1. Data loading)
data_err = load('task_2_errors.txt');

% 2. Columns selection
dt        = data_err(:, 1); 
err_euler = data_err(:, 2);
err_cn    = data_err(:, 3);

% 3. Log-Log Plot
figure('Name', 'Task 3 - Error Analysis sp', 'Color', 'w');


loglog(dt, err_euler, 'b-*', 'LineWidth', 1, 'MarkerSize', 8, ...
       'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b'); hold on;

loglog(dt, err_cn, 'r-*', 'LineWidth', 1, 'MarkerSize', 8, ...
       'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'r');

% 4. Labels
xlabel('Time Step \Delta t [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Absolute Error ||e_N||', 'FontSize', 12, 'FontWeight', 'bold');


% 5. Legend
legend('Implicit Euler', 'Crank-Nicolson', ...
       'Location', 'best', 'FontSize', 11);


grid on;
grid minor;
box on;