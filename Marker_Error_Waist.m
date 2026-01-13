function error = Marker_Error_Waist(params, rw_Waist_VIVE, rx_Waist_VIVE, ry_Waist_VIVE, rz_Waist_VIVE, pos_Waist_Original, time_column, z_Waist, x_Waist, y_Waist, t_corr)
% This function calculates the error between HVUT and VICON with physical constraints
% Integrates the original Marker_Error function with constraint handling
% 
% Parameters:
% params - optimization parameters [local_offset(1:3), distance_from_top_edge_to_base, width_of_UT_case, length_of_UT_case]
% rw, rx, ry, rz - rotation quaternion components
% pos_Original - original position data
% time_column - time data
% z_Waist, x_Waist, y_Waist - VICON position data
% t_corr - time correction offset

% Extract parameters
local_offset = params(1:3);
% distance_from_top_edge_to_base = params(4);
% width_of_UT_case = params(5);
% length_of_UT_case = params(6);

% Add 10mm margin for optimization flexibility
% h = distance_from_top_edge_to_base;
% w = width_of_UT_case;
% l = length_of_UT_case;

% Check if constraints are satisfied based on physical properties
% This could be implemented as nonlinear constraints in fmincon
% But here we add a penalty to the error if constraints are violated
% penalty = 0;
% if local_offset(1) < -l/2 || local_offset(1) > l/2
%     penalty = penalty + 1000;
% end
% if local_offset(2) < 0 || local_offset(2) > h
%     penalty = penalty + 1000;
% end
% if local_offset(3) < -w/2 || local_offset(3) > w/2
%    penalty = penalty + 1000;
% end

% Step 1: Add local offset and apply rotation
pos_Transformed = zeros(length(rw_Waist_VIVE), 3);
for i = 1:length(rw_Waist_VIVE)
    % Waist
    q = [rw_Waist_VIVE(i), rx_Waist_VIVE(i), ry_Waist_VIVE(i), rz_Waist_VIVE(i)];
    R = quat2rotm(q);
    pos_Transformed(i,:) = (R * (pos_Waist_Original(i,:) + local_offset)')';
end

% Step 2: Interpolation by Second
% Group data by second (based on time_column)
unique_seconds = unique(floor(time_column));
unique_seconds = unique_seconds(2:end-1); % Exclude first and last second

% Create containers for interpolated data
all_data_x = zeros(length(unique_seconds)*100, 1);
all_data_y = zeros(length(unique_seconds)*100, 1);
all_data_z = zeros(length(unique_seconds)*100, 1);

% Time frames array for plotting
time_frames = zeros(length(unique_seconds)*100, 1);
original_frame_counts = zeros(length(unique_seconds), 1);

% Process each second of data
for i = 1:length(unique_seconds)
    current_second = unique_seconds(i);
    indices = find(floor(time_column) == current_second);
    
    % Store original frame count for this second
    original_frame_counts(i) = length(indices);
    
    if length(indices) < 2
        continue; % 跳过这一秒的处理
    end
    
    % Extract the current second's data from transformed positions
    current_x_waist = pos_Transformed(indices, 1);
    current_y_waist = pos_Transformed(indices, 2);
    current_z_waist = pos_Transformed(indices, 3);
    
    % Create a normalized time vector for the current second (0 to 1)
    t_original = linspace(0, 1, length(indices));
    
    % Create a time vector for the desired 100 frames per second
    t_interp = linspace(0, 1, 100);
    
    % Interpolate data to have exactly 100 frames per second
    interp_x_waist = interp1(t_original, current_x_waist, t_interp, 'pchip');
    interp_y_waist = interp1(t_original, current_y_waist, t_interp, 'pchip');
    interp_z_waist = interp1(t_original, current_z_waist, t_interp, 'pchip');
    
    % Store interpolated data
    add_index = (i-1)*100+1 : i*100;
    all_data_x(add_index,1) = interp_x_waist';
    all_data_y(add_index,1) = interp_y_waist';
    all_data_z(add_index,1) = interp_z_waist';
    
    % Store the time frame
    time_frames(add_index,1) = ones(100, 1) * current_second;
end

% Step 3: Apply time offset correction (from get_time_offset)
% Define which VIVE points to compare based on t_corr
indexVICON = 600:min(length(z_Waist), length(all_data_y));
indexVIVE = indexVICON - t_corr;

% Make sure indices are within valid range
validIndices = (indexVIVE > 0) & (indexVIVE <= length(all_data_y)) & ...
               (indexVICON > 0) & (indexVICON <= length(z_Waist));
indexVICON = indexVICON(validIndices);
indexVIVE = indexVIVE(validIndices);

% Step 4: Apply Kabsch algorithm for best alignment
% Extract the corresponding points
VIVE_points = [-all_data_x(indexVIVE), all_data_z(indexVIVE), all_data_y(indexVIVE)];
VICON_points = [x_Waist(indexVICON), y_Waist(indexVICON), z_Waist(indexVICON)];

% Center the data points
VIVE_centroid = mean(VIVE_points, 1);
VICON_centroid = mean(VICON_points, 1);
VIVE_centered = VIVE_points - VIVE_centroid;
VICON_centered = VICON_points - VICON_centroid;

% Calculate the cross-covariance matrix
H = VIVE_centered' * VICON_centered;

% Calculate SVD of the cross-covariance matrix
[U, ~, V] = svd(H);

% Determine the rotation matrix
R_kabsch = V * U';

% Check for reflection
if det(R_kabsch) < 0
    V(:, 3) = -V(:, 3);
    R_kabsch = V * U';
end

% Calculate the translation vector
t_kabsch = VICON_centroid' - R_kabsch * VIVE_centroid';

% Transform VIVE points using the Kabsch transformation
VIVE_transformed = (R_kabsch * VIVE_points' + t_kabsch)';

% Step 5: Calculate RMS error between aligned VIVE and VICON data
errors_x = VIVE_transformed(:, 1) - VICON_points(:, 1);
errors_y = VIVE_transformed(:, 2) - VICON_points(:, 2);
errors_z = VIVE_transformed(:, 3) - VICON_points(:, 3);

% Total RMS error (combined across all dimensions)
rmsError = sqrt(mean(errors_x.^2 + errors_y.^2 + errors_z.^2));

% Return the error with penalty
error = rmsError;
end