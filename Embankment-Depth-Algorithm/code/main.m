%given csv file location here
csv_file = csvread('a240402b.csv');
ground_truth = readtable('ground_truth.xlsx');
x_length = 4;
y_length = 2;
%ground_truth = '';

%%% set the default values for peaks2 algorithm
minpeakheight = 0.1;
peaks2threshold = 0.001;
minpeakdistance = 0.02;

% Load CSV file
data_csv = csv_file; % Adjust the file name as needed

% Extract X, Y, and Z columns
X = data_csv(:, 1);
Y = data_csv(:, 2);
Z = data_csv(:, 3);

% Calculate the number of unique values in X and Y columns
num_unique_X = numel(unique(X));
num_unique_Y = numel(unique(Y));

% Find unique values in X and Y
unique_X = unique(X);
unique_Y = unique(Y);

% Find the second number in X and Y
second_X = unique_X(2);
second_Y = unique_Y(2);

% Convert to matrix format
[x, y] = meshgrid(linspace(0, x_length, num_unique_X), linspace(0, y_length, num_unique_Y));
z = griddata(X, Y, Z, x, y);

% Find peaks in the data
[pks, locs_y, locs_x] = peaks2(z, 'MinPeakHeight', minpeakheight, 'Threshold', peaks2threshold, 'MinPeakDistance',minpeakdistance);
%[pks, locs_y, locs_x] = peaks2(z, 'MinPeakHeight', 0.1, 'Threshold', 0.0015, 'MinPeakDistance',0.008);

% Plot the peak values in an X-Y grid using the scaled X values
scatter(second_X*locs_x, second_Y*locs_y, 80, pks, 'filled');
hold on; % Keep the current plot and add to it

% Read the groud truth Excel file
data_excel = ground_truth;

if isempty(ground_truth)
    disp('ground_truth is empty or null.');
    %saveas(gcf, '/Users/mohitsarin/Desktop/single_plot.png');
else
    % Extracting X, Y, and Z data
    x_excel = table2array(data_excel(:, 1));
    y_excel = table2array(data_excel(:, 2));
    z_excel = table2array(data_excel(:, 3));

    % Plot the scatter plot with 'X' markers
    scatter(x_excel, y_excel, 100, 'k', 'x'); % Increase marker size and use Z for color

    % Customize the plot
    xlabel('X-axis Label');
    ylabel('Y-axis Label');
    title('Combined Scatter Plot');
    colorbar; % Add color bar for the first scatter plot

    % Save the plot as an image
    saveas(gcf, '/Users/mohitsarin/Desktop/combined_plot.png');
end


%%%% error function is defined here: 
final_x_loc = second_X*locs_x;
final_y_loc = second_Y*locs_y;

% Convert to column vectors
final_x_loc = final_x_loc(:);
final_y_loc = final_y_loc(:);
final_z_loc = pks(:);
if isempty(ground_truth) 
    disp('ground_truth is empty or null.');
else
    x_excel = x_excel(:);
    y_excel = y_excel(:);
    z_excel = z_excel(:);
end

% Create matrices for predicted points and ground truth
predicted_points = [final_x_loc, final_y_loc];
if isempty(ground_truth)
    disp('ground_truth is empty or null.');
else
ground_truth = [x_excel, y_excel];
end

% Include all points from the input grid (used to find the peaks)
grid_points = [X, Y];

% Initialize arrays
assignments = zeros(size(ground_truth, 1), 1);
assigned_pred = zeros(size(predicted_points, 1), 1);

% Set threshold for assigning ground truth to predicted points
threshold = 0.3;
errors = zeros(size(ground_truth, 1), 1);
threshold_radius = 1.5

%addpath('/Users/mohitsarin/Desktop/TTI/finalcode/valley.m');
% Call the function to find valleys for each predicted point
valley_indices = find_valleys(predicted_points, grid_points, threshold);

% Print the predicted points and their corresponding valley points
for i = 1:size(predicted_points, 1)
    if valley_indices(i) ~= -1
        disp(['Predicted point: ', num2str(predicted_points(i, :)), ...
              ' -> Valley point: ', num2str(grid_points(valley_indices(i), :)), ...
              ' -> Valley height: ', num2str(Z(valley_indices(i)))]);
    else
        disp(['Predicted point: ', num2str(predicted_points(i, :)), ' -> No valley point found']);
    end
end

% Initialize array for valley heights
valley_heights = zeros(size(predicted_points, 1), 1);

% Calculate valley heights
for i = 1:size(predicted_points, 1)
    if valley_indices(i) ~= -1
        valley_heights(i) = Z(valley_indices(i)); % Assign Z value of valley point as valley height
    end
end

% Initialize arrays for embankment depth and percentage
embankment_depth = zeros(size(predicted_points, 1), 1);
embankment_percentage = zeros(size(predicted_points, 1), 1);
mat_thickness = 0.5; % Thickness of the material (assuming constant for all points)


% Print header
fprintf('%-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n', 'X_Predicted', 'Y_Predicted', 'X_Valley', 'Y_Valley', 'Max_Height', 'Valley_Height', 'Embankment_Depth', 'Embankment_Percentage');

% Calculate embankment depth and percentage for each predicted point
for i = 1:size(predicted_points, 1)
    if valley_indices(i) ~= -1
        max_height = final_z_loc(i); % Height at the predicted point
        valley_height = Z(valley_indices(i)); % Height at the valley point
        embankment_depth(i) = mat_thickness - (max_height - valley_height); % Calculate embankment depth
        embankment_percentage(i) = embankment_depth(i) / mat_thickness * 100; % Calculate embankment percentage
        fprintf('%-15.6f %-15.6f %-15.6f %-15.6f %-15.6f %-15.6f %-15.6f %-15.6f\n', ...
            predicted_points(i, 1), predicted_points(i, 2), ...
            grid_points(valley_indices(i), 1), grid_points(valley_indices(i), 2), ...
            max_height, valley_height, embankment_depth(i), embankment_percentage(i));
    else
        fprintf('%-15.6f %-15.6f %-15s %-15s %-15.6f %-15s %-15s %-15s\n', ...
            predicted_points(i, 1), predicted_points(i, 2), ...
            'N/A', 'N/A', max_height, 'N/A', 'N/A', 'N/A');
    end
end

if isempty(ground_truth)
    disp('ground_truth is empty or null.');
else

% Loop through each ground truth point
for i = 1:size(ground_truth, 1)
    % Calculate the Euclidean distances between the ground truth point and all unassigned predicted points
    dist_to_pred = sqrt(sum((predicted_points - ground_truth(i, :)).^2, 2));
    % Find the indices of nearby predicted points within the threshold
    nearby_indices = find(dist_to_pred <= threshold);

    if isempty(nearby_indices)
        assignments(i) = -1; % Mark as unassigned
        errors(i) = NaN;
    else
        % Find the index of the nearest predicted point that is not already assigned
        unassigned_indices = nearby_indices(assigned_pred(nearby_indices) == 0);
        if isempty(unassigned_indices)
            assignments(i) = -1; % Mark as unassigned if all nearby predicted points are already assigned
            errors(i) = NaN;
        else
            [~, nearest_idx] = min(dist_to_pred(unassigned_indices));
            nearest_idx = unassigned_indices(nearest_idx);
            assignments(i) = nearest_idx;
            assigned_pred(nearest_idx) = 1; 
            % Calculate error 
            errors(i) = sqrt(sum((predicted_points(nearest_idx, :) - ground_truth(i, :)).^2));
        end
    end
end

% Print the assignments with coordinates and errors
total_error = sum(errors(~isnan(errors)));
average_error = total_error / sum(~isnan(errors));
disp(['Total Euclidean error: ', num2str(total_error)]);
disp(['Average Euclidean error: ', num2str(average_error)]);

end