% =========================================================================
% Script MATLAB per Esercizio 2 - Convezione-Diffusione (Task 1A)
% =========================================================================
clear; clc; close all;

% 1. Imposta il nome del file da leggere
filename = 'task1a.txt';

% 2. Leggi i dati (saltando la prima riga di intestazione)
% readmatrix ignora automaticamente la prima riga di testo
data = readmatrix(filename, 'NumHeaderLines', 1);

% Rimuovi eventuali righe vuote (NaN) lasciate da Fortran per separare i blocchi
data = data(~any(isnan(data), 2), :);

% 3. Estrai le colonne nelle rispettive variabili
x       = data(:, 1);
t       = data(:, 2);
phi_num = data(:, 3);
phi_ex  = data(:, 4);

% 4. Trova tutti gli istanti di tempo unici salvati nel file
t_unique = unique(t);
num_times = length(t_unique);

% 5. Prepara la figura
figure('Name', 'Task 1A: Pure Convection', 'Position', [100, 100, 800, 500]);
hold on; grid on; box on;

% Crea una palette di colori per distinguere i vari tempi
colors = lines(num_times); 

% 6. Ciclo per plottare le curve ad ogni istante di tempo
for i = 1:num_times
    % Trova gli indici dei dati corrispondenti al tempo attuale
    idx = (t == t_unique(i));
    
    % Estrai i segmenti per questo istante
    x_plot       = x(idx);
    phi_num_plot = phi_num(idx);
    phi_ex_plot  = phi_ex(idx);
    
    % Plot Soluzione Esatta (Linea continua spessa)
    plot(x_plot, phi_ex_plot, '-', 'Color', colors(i,:), 'LineWidth', 1, ...
        'DisplayName', sprintf('Exact (t = %.2f)', t_unique(i)));
    
    % Plot Soluzione Numerica (Pallini colorati)
    plot(x_plot, phi_num_plot, 'o', 'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), 'MarkerSize', 2, ...
        'DisplayName', sprintf('Num CS2 (t = %.2f)', t_unique(i)));
end

% 7. Formattazione estetica del grafico
title('Task 1A: Pure Convection - CS2 Scheme (\nu = 0, c = 0.1)', 'FontSize', 12);
xlabel('Spatial Coordinate, x', 'FontSize', 11);
ylabel('\phi(x,t)', 'FontSize', 11);

% Limiti degli assi per inquadrare bene la Gaussiana
xlim([0, 1]);
ylim([-0.1, 1.2]); 

% Posiziona la legenda fuori dal grafico a destra per non coprire l'onda
legend('Location', 'eastoutside', 'FontSize', 9);

hold off;