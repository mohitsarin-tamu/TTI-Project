% Load CSV file
data = csvread('/Users/mohitsarin/Desktop/a240402b.csv'); % Adjust the file name as needed

% Extract X, Y, and Z columns
X = data(:, 1);
Y = data(:, 2);
Z = data(:, 3);

% Convert to matrix format
[x, y] = meshgrid(linspace(min(X), max(X), (max(X)- min(X))/0.011811000000000), linspace(min(Y), max(Y), (max(Y) -min(Y))/0.019444000000000));
z = griddata(X, Y, Z, x, y);

% Save the griddata result in matrix format
save('/Users/mohitsarin/Desktop/griddata_result.mat', 'z');
