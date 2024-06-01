% Load CSV file
data = csvread('/Users/mohitsarin/Desktop/a240402b.csv'); % Adjust the file name as needed

% Extract X, Y, and Z columns
X = data(:, 1);
Y = data(:, 2);
Z = data(:, 3);

% Convert to matrix format
[x, y] = meshgrid(linspace(min(X), max(X), 10), linspace(min(Y), max(Y), 10));
z = griddata(X, Y, Z, x, y);

% % Find peaks
% ix = find(imregionalmax(z, 8));

% % Plot surface
% surf(x, y, z, 'FaceColor', 'interp');
% hold on;
% 
% % Plot peaks
% plot3(x(ix), y(ix), z(ix), 'r*', 'MarkerSize', 24);

% Save the image in the same location
saveas(gcf, '/Users/mohitsarin/Desktop/plot.png'); % Adjust the file name as needed
