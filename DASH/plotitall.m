% overlay_pressure_recordings.m
% Overlays two specific pressure recording CSV files on the same plot.
%
% Assumptions:
%   Column 1 = time (ms)
%   Column 2 = pressure (psi)

clc; clear; close all;

% -------- FILES TO PLOT --------
fileList = {
    "Pressure Recording 2025-09-01_173604.csv"
    "Pressure Recording 2025-09-01_172352.csv"
};

colors = {'b','r'}; % Blue and red for clarity

figure('Name','Overlayed Pressure Recordings','NumberTitle','off'); 
hold on;

for k = 1:numel(fileList)
    % Read CSV file
    data = readmatrix(fileList{k});
    
    % Extract time and pressure
    time_ms = data(:,1);
    pressure_psi = data(:,2);
    
    % Plot
    plot(time_ms/1000, pressure_psi, 'Color', colors{k}, 'LineWidth', 1.5, ...
         'DisplayName', strrep(fileList{k}, '_', '\_'));
end

% Labels and legend
xlabel('Time (s)');
ylabel('Pressure (psi)');
title('Extended Duration Test');
grid on;
legend('Location','best');
