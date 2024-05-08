% Read the Excel file
[~, ~, data] = xlsread('ground_truth.xlsx');

% Convert cell array to numeric array
data = cell2mat(data(2:end, :));

% Extracting X and Y data
x = data(:, 1);
y = data(:, 2);

% Plot the scatter plot
scatter(x, y, 50, 'filled',"+");

% Customize the plot
xlabel('X-axis Label');
ylabel('Y-axis Label');
title('Scatter Plot');
ylim([0, 0.5]); % Set the y-axis limit to 0.5

% Display the plot
fig = gcf; % Get current figure handle

% Save the plot as an image
saveas(fig, 'scatter_plot.png');

% Close the figure window
% close(fig);
