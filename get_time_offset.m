function t_corr = get_time_offset(z_Waist, all_data_y_waist, z_RToe, all_data_y_rtoe, z_LToe, all_data_y_ltoe, graph)

% Sript to find time synchronisation

% VICON data to use
indexVICON = 600:min(length(z_Waist), length(all_data_y_waist));
error = zeros(1,599); 

for t_offset = 1:599
    % define which HVUT points to compare
    indexHVUT = indexVICON - t_offset; 

    % calculate average vertical offset between systems
    global_offset = mean(z_Waist(indexVICON)) - mean(all_data_y_waist(indexHVUT)); 

    % Get RMS error between VIVE and VICON
    error(t_offset) = sqrt(sum((z_Waist(indexVICON) - all_data_y_waist(indexHVUT) - global_offset).^2));
end

% use minimum RMS error to calculate time correction
[~, t_corr] = min(error);

% if graph
%     indexHVUT = indexVICON - t_corr;
%     global_offset = mean(z_Waist(indexVICON)) - mean(all_data_y_waist(indexHVUT));
%     
%     figure
%     hold on
%     plot(indexVICON, z_Waist(indexVICON),'b-')
%     indexHVUT = indexVICON - t_corr;
%     plot(indexVICON, all_data_y_waist(indexHVUT) + global_offset, 'r-')
% end

if graph
    indexHVUT = indexVICON - t_corr;
    
    % Waist comparison (original)
    global_offset_waist = mean(z_Waist(indexVICON)) - mean(all_data_y_waist(indexHVUT));
    figure
    hold on
    plot(indexVICON, z_Waist(indexVICON),'b-', 'DisplayName', 'Vicon Waist Z')
    plot(indexVICON, all_data_y_waist(indexHVUT) + global_offset_waist, 'r-', 'DisplayName', 'VIVE Waist Y + offset')
    title('Waist Synchronization (Vicon Z vs VIVE Y)')
    xlabel('Frame Index')
    ylabel('Position (mm)')
    legend
    grid on
    
    % RToe comparison
    global_offset_rtoe = mean(z_RToe(indexVICON)) - mean(all_data_y_rtoe(indexHVUT));
    figure
    hold on
    plot(indexVICON, z_RToe(indexVICON),'b-', 'DisplayName', 'Vicon RToe Z')
    plot(indexVICON, all_data_y_rtoe(indexHVUT) + global_offset_rtoe, 'r-', 'DisplayName', 'VIVE RToe Y + offset')
    title('Right Toe Synchronization (Vicon Z vs VIVE Y)')
    xlabel('Frame Index')
    ylabel('Position (mm)')
    legend
    grid on
    
    % LToe comparison
    global_offset_ltoe = mean(z_LToe(indexVICON)) - mean(all_data_y_ltoe(indexHVUT));
    figure
    hold on
    plot(indexVICON, z_LToe(indexVICON),'b-', 'DisplayName', 'Vicon LToe Z')
    plot(indexVICON, all_data_y_ltoe(indexHVUT) + global_offset_ltoe, 'r-', 'DisplayName', 'VIVE LToe Y + offset')
    title('Left Toe Synchronization (Vicon Z vs VIVE Y)')
    xlabel('Frame Index')
    ylabel('Position (mm)')
    legend
    grid on
end