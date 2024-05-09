% Read the Excel file
data = readtable('/Users/mohitsarin/Desktop/TTI/ground_truth.xlsx')
% Convert table to array
data = table2array(data);

% Extracting X, Y, and Z data
x = data(:, 1);
y = data(:, 2);
z = data(:, 3);

% Plot the scatter plot
%scatter(x, y, 100, z, 'filled'); % Increase marker size and use Z for color


% Plot the scatter plot with 'X' markers
scatter(x, y, 100,'x'); % Increase marker size and use Z for color

% Customize the plot
xlabel('X-axis Label');
ylabel('Y-axis Label');
title('Scatter Plot');
colorbar; % Add color bar for Z values

% Save the plot as an image
saveas(gcf, '/Users/mohitsarin/Desktop/base_plot.png');
