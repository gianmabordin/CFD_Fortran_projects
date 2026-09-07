% =========================================================================
% Script MATLAB per Task 2 - Pure Diffusion (CS4)
% =========================================================================
clear; clc; close all;

filename = 'task2.txt';

% 1. Caricamento dati
if ~exist(filename, 'file')
    error('File %s non trovato. Lancia prima il codice Fortran!', filename);
end

data = readmatrix(filename, 'NumHeaderLines', 1);
data = data(~any(isnan(data), 2), :); % Rimuove eventuali righe vuote o NaN

% 2. Estrazione colonne
x       = data(:, 1);
t       = data(:, 2);
phi_num = data(:, 3);
phi_ex  = data(:, 4);

t_unique = unique(t);
num_times = length(t_unique);

% 3. Creazione Figura
figure('Name', 'Task 2: Pure Diffusion - CS4', 'Color', 'w');
hold on; grid on; box on;

% Palette colori sfumata per mostrare il passare del tempo
colors = winter(num_times); 

for i = 1:num_times
    idx = (t == t_unique(i));
    
    % Soluzione Esatta (Linea continua)
    plot(x(idx), phi_ex(idx), '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Exact (t = %.1f)', t_unique(i)));
    
    % Soluzione Numerica (Cerchietti)
    plot(x(idx), phi_num(idx), 'o', 'MarkerEdgeColor', colors(i,:), ...
        'MarkerSize', 4, 'DisplayName', sprintf('CS4 (t = %.1f)', t_unique(i)));
end

% 4. Formattazione Estetica
title(['Task 2: Pure Diffusion (\nu = 0.01) - CS4 Scheme'], 'FontSize', 12);
xlabel('Spatial Coordinate x', 'FontSize', 11);
ylabel('\phi(x,t)', 'FontSize', 11);
xlim([0, 1]);
ylim([-0.1, 1.1]); % Manteniamo il range visibile
legend('Location', 'northeastoutside');

hold off;