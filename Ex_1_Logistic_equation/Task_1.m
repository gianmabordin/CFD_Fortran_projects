% ==========================================
% PLOT TASK 1: METHOD COMPARISON FOR DT = 1
% ==========================================
clear; clc; close all;

% 1. Data loading
dati_euler = load('task_1eul.txt');
dati_cn    = load('task_1cn.txt');

% 2. Column selection
t         = dati_euler(:, 1); 
phi_exact = dati_euler(:, 3); 
phi_euler = dati_euler(:, 2);
phi_cn    = dati_cn(:, 2);

% 3. Plot
figure('Name', 'Task 1', 'Color', 'w');

% Curves
plot(t, phi_euler, 'b-*', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on; 
plot(t, phi_cn, 'r-*', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');             
plot(t, phi_exact, 'k--', 'LineWidth', 1, 'MarkerSize', 4, 'MarkerEdgeColor', 'k');         

% Labels
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Solution \phi(t)', 'FontSize', 12, 'FontWeight', 'bold');

% Legend
legend('Implicit Euler', 'Crank-Nicolson', 'Analytical solution', 'Location', 'best', 'FontSize', 11);

grid on;      
grid minor;  
box on;       