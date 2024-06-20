% Function to find valleys for each predicted point
function [valley_indices] = find_valleys(predicted_points, grid_points, threshold_radius)
    % Initialize the output variable
    valley_indices = zeros(size(predicted_points, 1), 1);

    % Loop through each predicted point
    for i = 1:size(predicted_points, 1)
        % Calculate the Euclidean distances between the current predicted point and all grid points
        dist_to_all_points = sqrt(sum((grid_points - predicted_points(i, :)).^2, 2));
        
        % Find the indices of nearby points within the search radius
        nearby_indices = find(dist_to_all_points <= threshold_radius);
        
        % Sort the distances to get the nearest points
        [~, sorted_indices] = sort(dist_to_all_points(nearby_indices));
        
        % Extract the distances of nearby points
        sorted_distances = dist_to_all_points(nearby_indices(sorted_indices));
        
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
