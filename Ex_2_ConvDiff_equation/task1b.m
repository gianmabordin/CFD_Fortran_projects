% =========================================================================
% Script MATLAB per Esercizio 2 - Convezione-Diffusione (Task 1B)
% =========================================================================
clear; clc; close all;

% Definisci i file generati da Fortran e i nomi degli schemi
filenames = {'task1b_fw1.txt', 'task1b_bw1.txt', 'task1b_cs2.txt'};
schemes   = {'FW1 (Forward 1st Order)', 'BW1 (Backward 1st Order)', 'CS2 (Centered 2nd Order)'};

% Ciclo automatico per creare i 3 grafici
for k = 1:length(filenames)
    
    % 1. Leggi i dati (saltando l'intestazione e togliendo i NaN)
    data = readmatrix(filenames{k}, 'NumHeaderLines', 1);
    data = data(~any(isnan(data), 2), :);
    
    % 2. Estrai le colonne
    x       = data(:, 1);
    t       = data(:, 2);
    phi_num = data(:, 3);
    phi_ex  = data(:, 4);
    
    % 3. Trova i tempi unici salvati
    t_unique = unique(t);
    num_times = length(t_unique);
    
    % 4. Crea una nuova figura per ogni schema
    figure('Name', sprintf('Task 1B: %s', schemes{k}), 'Position', [100+(k*50), 100+(k*50), 800, 500]);
    hold on; grid on; box on;
    
    colors = lines(num_times); 
    
    % 5. Plotta le curve per ogni istante di tempo
    for i = 1:num_times
        idx = (t == t_unique(i));
        
        x_plot       = x(idx);
        phi_num_plot = phi_num(idx);
        phi_ex_plot  = phi_ex(idx);
        
        % Soluzione Esatta (Linea continua)
        plot(x_plot, phi_ex_plot, '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Exact (t = %.0f)', t_unique(i)));
        
        % Soluzione Numerica (Pallini)
        plot(x_plot, phi_num_plot, 'o', 'Color', colors(i,:), ...
            'MarkerFaceColor', colors(i,:), 'MarkerSize', 4, ...
            'DisplayName', sprintf('Num (t = %.0f)', t_unique(i)));
    end
    
    % 6. Formattazione specifica per l'onda Armonica
    title(sprintf('Task 1B: Pure Convection - %s', schemes{k}), 'FontSize', 12);
    xlabel('Spatial Coordinate, x', 'FontSize', 11);
    ylabel('\phi(x,t)', 'FontSize', 11);
    
    % Limiti asse Y adattati per l'onda sinusoidale (da -1 a 1)
    xlim([0, 1]);
    ylim([-1.5, 1.5]); 
    
    legend('Location', 'eastoutside', 'FontSize', 9);
    hold off;
    
end