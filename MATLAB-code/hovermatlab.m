% Load CSV file
data_csv = csvread('/Users/mohitsarin/Desktop/TTI/a240402b.csv'); % Adjust the file name as needed

% Extract X, Y, and Z columns
X = data_csv(:, 1);
Y = data_csv(:, 2);
Z = data_csv(:, 3);

figure;
scatter3(X, Y, Z, 'b.');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Scatter Plot with Hover Information');
grid on;

% Enable hovering
dcm_obj = datacursormode(gcf);
set(dcm_obj, 'UpdateFcn', @hover_callback);

function txt = hover_callback(~, event_obj)
    pos = get(event_obj, 'Position');
    txt = {['X: ', num2str(pos(1))], ['Y: ', num2str(pos(2))], ['Z: ', num2str(pos(3))]};
end
