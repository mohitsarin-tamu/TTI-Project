function [pks,locs_y,locs_x]=peaks2(data,varargin)
% Find local peaks in 2D data.
% Syntax chosen to be as close as possible to the original Matlab
% 'findpeaks' function but not require any additional toolbox.
%
% SYNTAX:
% pks=peaks(data) finds local peaks.
%
% [pks,locs_y,locs_x]=peaks(data) finds local peaks and their array
% coordinates.
%
% [pks,locs_y,locs_x]=peaks(...,'MinPeakHeight',{scalar value}) only retains
% those peaks which are equal to or greater than this absolute value.
%
% [pks,locs_y,locs_x]=peaks(...,'Threshold',{scalar value}) only retains
% those peaks that are higher than their immediate surroundings by this value.
%
% [pks,locs_y,locs_x]=peaks(...,'MinPeakDistance',{scalar value}) finds peaks
% separated by more than the specified minimum CARTESIAN peak distance (a
% circle around the peak). It starts from the strongest peak and goes
% iteratively lower. Any peak 'shadowed' in the vicinity of a stronger
% peak is discarded.
%
% ALGORITHM: 
% A peak is considered to be a data point strictly greater than its
% immediate neighbors. You can change this condition to 'greater or equal'
% in the code, but be aware that in such case, it might create false
% detections in flat areas, but these can be accounted for by introduction
% of a small Threshold value.
%
% Even though this function is shared here free for use and any
% modifications you might find useful, I would appreciate if you would
% quote me in case you are going to use this function for any non-personal
% tasks.
% (C) Kristupas Tikuisis 2023.

% data = '/Users/mohitsarin/Desktop/a240402b_data.mat'
%% Initial data check
% Let's simplify the function. Let it work on 1D or 2D data only, and check
% if the input data satisfies this criteria:
if ~ismatrix(data)
    error('Only 1D (vectors) or 2D (matrices) data accepted.');
end


%% Locate all peaks
% A peak is a data point HIGHER than its immediate neighbors. There are 8
% around eaxh point, and we will go through each. Oh yes, no escaping that.
%
% To introduce as little of intermediate variables and keep their memory
% footprint as low as possible, I will introduce 2 logical arrays: one to
% mark the peaks and be iteratively updated until we check all its
% neighbors; and another temporal variable just to prepare data for
% comparison (mind about the edge points which do not have any neighbors!):
ispeak=false(size(data)); % to store peak flags.
isgreater=true(size(data)); % to store comparison result for one particular neighbor.

% Now start analyzing every data point.
%
% 1st neighbor immediatelly to the left:
ispeak=([true(size(data,1),1) [data(:,2:end)>data(:,1:end-1)]]); % for this case, we can update the peak array directly.
%
% 2nd neighbor at the top-left:
isgreater(2:end,2:end)=(data(2:end,2:end)>data(1:end-1,1:end-1)); % this time, due to points on the diagonal, a temporary array will have to be involved.
ispeak=ispeak&isgreater; % time to update the peak array.
%
% 3rd neighbor immediatelly at the top:
ispeak=ispeak&([true(1,size(data,2)); (data(2:end,:)>data(1:end-1,:))]); % once again, for this case, we can update the peak array directly.
%
% 4th neighbor at the top-right:
isgreater=true(size(data)); % rebuild a fresh array for comparison...
isgreater(2:end,1:end-1)=(data(2:end,1:end-1)>data(1:end-1,2:end));
ispeak=ispeak&isgreater;
%
% 5th neighbor immediatelly to the right:
ispeak=ispeak&([(data(:,1:end-1)>data(:,2:end)) true(size(data,1),1)]); % once again, for this case, we can update the peak array directly.
%
% 6th neighbor to the bottom-right:
isgreater=true(size(data));
isgreater(1:end-1,1:end-1)=(data(1:end-1,1:end-1)>data(2:end,2:end));
ispeak=ispeak&isgreater;
%
% 7th neighbor immediatelly at the bottom:
ispeak=ispeak&([(data(1:end-1,:)>data(2:end,:)); true(1,size(data,2))]); % once again, for this case, we can update the peak array directly.
%
% 8th neighbor at the bottom-left:
isgreater=true(size(data));
isgreater(1:end-1,2:end)=(data(1:end-1,2:end)>data(2:end,1:end-1));
ispeak=ispeak&isgreater;
%
% Discard the temporary variable:
clear isgreater
% By now, raw peak indentification is completed.


%% Return final results
% First, essential raw peak location steps. Peak values and locations:
locs=find(ispeak); %clear ispeak
pks=data(locs);
% Note that LINEAR indices have been returned. Let's turn them to array
% indices for final output:
[locs_y,locs_x]=ind2sub(size(data),locs);


%% Perform customary post-processing determined by optional function parameters.
% First, check if any parameters were supplied:
if isempty(varargin)
    return
end
% ...then a quality check - the parameters should come in pairs, therefore
% the length should be even:
if mod(length(varargin),2)~=0
    warning('Optional name-value parameters should come in pairs. Something is missing. Peak search will proceed with default values.');
    return
end
% ...after this step, we can sort the input into names (parameters) and
% values:
params=varargin(1:2:end);
vals=varargin(2:2:end);
% ...last quality check - all even values should be CHAR entries specifying
% a parameter to be adjusted:
ischarparam=cellfun(@ischar,params);
if any(~ischarparam)
    warning('Some parameters are not written as characters and not recognisable. They should always come is name(char)-value pairs. Peak search will proceed with default values.');
    return
end; clear ischarparam varargin

%--------------------------------------------------------------------------
% No go over all supplied parameters and check the found peaks accordingly.

% 1. MinPeakHeight - absolute minimum value for a peak:
isparm=find(cellfun(@(x)isequal(x,'MinPeakHeight'),params),1);
if ~isempty(isparm)
    % Locate which peaks satisfy this condition:
    suitable=(pks>=vals{isparm});
    % ...and only keep those:
    pks=pks(suitable);
    locs=locs(suitable);
    clear suitable
end

% 2. Threshold - peak must be greater than its neighbours by this value.
isparm=find(cellfun(@(x)isequal(x,'Threshold'),params),1);
if ~isempty(isparm)
    % For this, we will need to convert linear indices to array indices:
    [row_y,col_x]=ind2sub(size(data),locs);
    % These will be the original indices.
    
    % This is how array coordinates would change relatively around each
    % peak (y,x):
    % (-1,-1)  (-1,0)  (-1,+1)
    % ( 0,-1)  ( 0,0)  ( 0,+1)
    % (+1,-1)  (+1,0)  (+1,+1)
    % ...turned into vectors disregarding the (0,0), the center data point:
    delta_y=[-1 -1 -1 0 0 +1 +1 +1];
    delta_x=[-1 0 +1 -1 +1 -1 0 +1];
    % Let's add these deltas to the detected peak positions to get the
    % coordinates of their immediate neighbors:
    neighbor_locs_y=row_y+delta_y;
    neighbor_locs_x=col_x+delta_x; clear row_y col_x delta_x delta_y
    % ...don't forget to check for unrealistic indices beyond array
    % borders:
    neighbor_locs_y(neighbor_locs_y<1)=1;
    neighbor_locs_y(neighbor_locs_y>size(data,1))=size(data,1);
    neighbor_locs_x(neighbor_locs_x<1)=1;
    neighbor_locs_x(neighbor_locs_x>size(data,2))=size(data,2);
    % ...convert to linear indices:
    neighbor_locs=sub2ind(size(data),neighbor_locs_y,neighbor_locs_x);
    clear neighbor_locs_y neighbor_locs_x
    
    % So we have neighbor values by now. Are our peaks higher than those by
    % the set Threshold value?
    suitable=(data(locs)-vals{isparm}>=data(neighbor_locs));
    % Now check for those cases when by mistake (earlier step for checking
    % for indices beyond array boundaries) a peak itself is taken as a
    % neighbor as well:
    suitable(data(neighbor_locs)==data(locs))=true;
    
    % Only keep those is they are greater than ALL neighbors (in other
    % words, those elements where NONE are lesser):
    suitable=~any(~suitable,2); clear neighbor_locs
    
    % Final step - locate suitable element indices:
    suitable=find(suitable);
    
    % That's it, modify the output array:
    pks=pks(suitable);
    locs=locs(suitable); clear suitable
end

% 3. 'MinPeakDistance'
isparm=find(cellfun(@(x)isequal(x,'MinPeakDistance'),params),1);
if ~isempty(isparm)
    
    % First, sort the peaks in order of amplitude:
    [pks_sorted,idx]=sort(pks,'descend');
    locs_sorted=locs(idx); clear idx

    % The flow is as follows: start from the highest peak and discard any
    % other peaks closer than the CARTESIAN MinPeakDistance (that is, the
    % CIRCLE around the peak is going to be checked); then continue until
    % the whole list (updated iteratively as items might get removed) has
    % been checked.
    
    % Convert locations to array indices:
    [row_y,col_x]=ind2sub(size(data),locs_sorted);

    % Start from the highest peak:
    this_peak=1;
    while this_peak<(length(pks_sorted)+1)

        % Cartesian distances to ALL its remaining & yet unchecked neighbors
        % (including itself):
        dist=sqrt((row_y-row_y(this_peak)).^2+(col_x-col_x(this_peak)).^2);
        % Now simply check which neighbors are WITHIN the MinPeakDistance
        % BUT also nonzero (the peak should not be compared to itself):
        within=( (dist<=vals{isparm}) & (dist~=0) );
        % ...and delete those entries satisfying the condition:
        pks_sorted(within)=[];
        locs_sorted(within)=[];
        row_y(within)=[]; col_x(within)=[];

        % Update the peak counter:
        this_peak=this_peak+1;

    end; clear this_peak within dist row_y col_x

    % Update the peak location and value list:
    pks=pks_sorted;
    locs=locs_sorted;
end

% I haven't figured out how to do it more elegantly, but turn the default
% linear indices to array indices:
[locs_y,locs_x]=ind2sub(size(data),locs);

%==========================================================================
% End of the entire function

end


% Load CSV file
data_csv = csvread('/Users/mohitsarin/Desktop/TTI/a240402b.csv'); % Adjust the file name as needed

% Extract X, Y, and Z columns
X = data_csv(:, 1);
Y = data_csv(:, 2);
Z = data_csv(:, 3);

% Scale X and Y (if needed)
X_scaled = X;
Y_scaled = Y;

% Convert to matrix format
[x, y] = meshgrid(linspace(0, 4, 338), linspace(0, 2, 104));
z = griddata(X, Y, Z, x, y);

% Find peaks in the data
[pks, locs_y, locs_x] = peaks2(z, 'MinPeakHeight', 0.1, 'Threshold', 0.001, 'MinPeakDistance',0.02);
%[pks, locs_y, locs_x] = peaks2(z, 'MinPeakHeight', 0.1, 'Threshold', 0.0015, 'MinPeakDistance',0.008);

% Plot the peak values in an X-Y grid using the scaled X values
scatter(0.0118*locs_x, 0.0194*locs_y, 80, pks, 'filled');
hold on; % Keep the current plot and add to it

% Read the Excel file
data_excel = readtable('/Users/mohitsarin/Desktop/TTI/ground_truth.xlsx');

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


%%%% error function is defined here: 
final_x_loc = 0.0118*locs_x;
final_y_loc = 0.0194*locs_y;

%% distance with ground truth:

%distance = (final_x_loc - x_excel)**2  + (final_y_loc - y_excel) **2;

%%%%%%
% Convert to column vectors
final_x_loc = final_x_loc(:);
final_y_loc = final_y_loc(:);
x_excel = x_excel(:);
y_excel = y_excel(:);
% new
z_excel = z_excel(:);
final_z_loc = pks(:);

% Create matrices for predicted points and ground truth
predicted_points = [final_x_loc, final_y_loc];
ground_truth = [x_excel, y_excel];
  
% Initialize arrays
assignments = zeros(size(ground_truth, 1), 1);
assigned_pred = zeros(size(predicted_points, 1), 1);

% Set threshold for assigning ground truth to predicted points
threshold = 0.3;

errors = zeros(size(ground_truth, 1), 1);

% valley points detection for each predicted point - check the valley
% points for the nearby peak

function [valley_indices] = find_valleys(predicted_points)
    % Initialize the output variable
    valley_indices = zeros(size(predicted_points, 1), 1);

    % Loop through each predicted point
    for i = 1:size(predicted_points, 1)
        % Calculate the Euclidean distances between the current predicted point and all other predicted points
        dist_to_other_pred = sqrt(sum((predicted_points - predicted_points(i, :)).^2, 2));
        
        % Find the indices of nearby predicted points within the search radius
        nearby_indices = find(dist_to_other_pred <= 1.0); % DEFINING the limit for valley indices
        
        % Exclude the current predicted point itself from nearby indices
        nearby_indices(nearby_indices == i) = [];

        % Sort the distances to get the nearest points
        [~, sorted_indices] = sort(dist_to_other_pred(nearby_indices));
        
        % Extract the distances of nearby points
        sorted_distances = dist_to_other_pred(nearby_indices(sorted_indices));
        
        % Calculate the 5th lowest percentile distance
        fifth_percentile_index = ceil(0.05 * numel(sorted_distances));
        threshold_distance = sorted_distances(fifth_percentile_index);
        
        % Find the indices of points that are within the threshold distance
        within_threshold_indices = nearby_indices(sorted_indices(sorted_distances <= threshold_distance));
        
        % If there are no points within the threshold, assign -1
        if isempty(within_threshold_indices)
            valley_indices(i) = -1;
        else
            % Find the index of the nearest point among the points within the threshold
            nearest_index = within_threshold_indices(1);
            valley_indices(i) = nearest_index;
        end
    end
end

% Call the function to find valleys for each predicted point
valley_indices = find_valleys(predicted_points);

% Loop through each ground truth point
for i = 1:size(ground_truth, 1)
    % Calculate the Euclidean distances between the ground truth point and all unassigned predicted points
    dist_to_pred = sqrt(sum((predicted_points - ground_truth(i, :)).^2, 2));
    % Find the indices of nearby predicted points within the threshold
    nearby_indices = find(dist_to_pred <= threshold);
    %disp(dist_to_pred);
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

for i = 1:length(assignments)
    if assignments(i) == -1
        disp(['Ground truth point ', num2str(i), ' (x:', num2str(ground_truth(i,1)), ', y:', num2str(ground_truth(i,2)), ') is unassigned']);
    else
        disp(['Ground truth point ', num2str(i), ' (x:', num2str(ground_truth(i,1)), ', y:', num2str(ground_truth(i,2)), ') assigned to predicted point ', num2str(assignments(i)), ' (x:', num2str(predicted_points(assignments(i),1)), ', y:', num2str(predicted_points(assignments(i),2)), '), error: ', num2str(errors(i))]);
        % Print valley point for this assigned predicted point
        assigned_valley_idx = valley_indices(assignments(i));
        if assigned_valley_idx ~= -1 % If a valley point is found
            disp(['  Valley point: (x:', num2str(predicted_points(assigned_valley_idx,1)), ', y:', num2str(predicted_points(assigned_valley_idx,2)), ')']);
        else
            disp('  No valley point found');
        end
    end
end

% Step by step code to show unassigned predicted points
disp('Unassigned Predicted Points:');
unassigned_indices = find(assigned_pred == 0);
for i = 1:length(unassigned_indices)
    disp(['Predicted point ', num2str(unassigned_indices(i)), ' (x:', num2str(predicted_points(unassigned_indices(i),1)), ', y:', num2str(predicted_points(unassigned_indices(i),2)), ') is unassigned']);
end

% Calculate embankment depth and embankment percentage
embankment_depth = zeros(length(assignments), 1);
embankment_percentage = zeros(length(assignments), 1);
mat_thickness = 0.5

for i = 1:length(assignments)
    if assignments(i) ~= -1 % If the ground truth point is assigned
        % Calculate embankment depth
        max_height = final_z_loc(assignments(i)); % Max height from predicted points
        valley_height = final_z_loc(valley_indices(assignments(i))); % Valley height from predicted points
        embankment_depth(i) = mat_thickness - (max_height - valley_height);
        % Calculate embankment percentage
        embankment_percentage(i) = embankment_depth(i) / mat_thickness * 100;
    end
end

% Display embankment depth and embankment percentage for assigned ground truth points
for i = 1:length(assignments)
    if assignments(i) ~= -1 % If the ground truth point is assigned
        disp(['Embankment depth for ground truth point ', num2str(i), ': ', num2str(embankment_depth(i))]);
        disp(['Embankment percentage for ground truth point ', num2str(i), ': ', num2str(embankment_percentage(i)), '%']);
    end
end

% Filter out unassigned ground truth points
assigned_indices = find(assignments ~= -1);
embankment_percentage_assigned = embankment_percentage(assigned_indices);

% Plot embankment percentage against assigned ground truth points
figure;
plot(assigned_indices, embankment_percentage_assigned, '-o', 'LineWidth', 2);
xlabel('Ground Truth Point Index');
ylabel('Embankment Percentage');
title('Embankment Percentage vs. Ground Truth Point Index');
grid on;




% Call the function to find valleys for each predicted point
% valley_indices = find_valleys(predicted_points);

% There are few points we have to find z values for this.
