%ALIGNDIMENSIONS create a new reference array to align a new grid to an existing one while preserving the target cell size

function [new_XWorldLimits, new_YWorldLimits, x_offset_cells, y_offset_cells, width_cells, height_cells] = alignDimensions(target_XWorldLimits, target_YWorldLimits, original_XWorldLimits, original_YWorldLimits, x_celldim, y_celldim, x_origin, y_origin)

%     original_x_min = XWorldLimits(1); %max DEM X
%     original_x_max = XWorldLimits(2); %min DEM X
%     original_y_min = YWorldLimits(1); %max DEM Y
%     original_y_max = YWorldLimits(2); %min DEM Y
% 
%     % Define the original grid's X and Y world limits
%     original_XWorldLimits = [original_x_min, original_x_max];
%     original_YWorldLimits = [original_y_min, original_y_max];

    
    % Calculate the dimensions of the input grid
    width_cells = round((target_XWorldLimits(2) - target_XWorldLimits(1)) / abs(x_celldim));
    height_cells = round((target_YWorldLimits(2) - target_YWorldLimits(1)) / abs(y_celldim));
    
    % Determine the alignment corner
    % x_origin = info.SpatialRef.RowsStartFrom;  % Set your desired alignment ("west" or "east")
    % y_origin = info.SpatialRef.ColumnsStartFrom; % Set your desired alignment ("north" or "south")
    
    % Calculate the offset in cells based on the alignment
    switch x_origin
        case "west"
            x_offset_cells = round((-original_XWorldLimits(1) + target_XWorldLimits(1)) / x_celldim);
            x_offset_cells = [x_offset_cells x_offset_cells+width_cells];
        case "east"
            x_offset_cells = round((original_XWorldLimits(2) - target_XWorldLimits(2)) / x_celldim) - width_cells; %not tested
            x_offset_cells = [x_offset_cells-width_cells x_offset_cells];
    end
    
    switch y_origin
        case "north"
            y_offset_cells = round((-original_YWorldLimits(2) + target_YWorldLimits(2)) / y_celldim);
            y_offset_cells = [y_offset_cells y_offset_cells+height_cells];
        case "south"
            y_offset_cells = round((original_YWorldLimits(2) - target_YWorldLimits(2)) / y_celldim) - height_cells; %not tested
            y_offset_cells = [y_offset_cells-height_cells y_offset_cells];
    end

    
    % Calculate the new world limits based on the specified origin and alignment
    new_x_min = original_XWorldLimits(1) + x_offset_cells(1) * x_celldim;
    new_x_max = original_XWorldLimits(1) + x_offset_cells(2) * x_celldim;
    new_y_min = original_YWorldLimits(2) + y_offset_cells(1) * y_celldim;
    new_y_max = original_YWorldLimits(2) + y_offset_cells(2) * y_celldim;

    new_XWorldLimits = [min(new_x_min, new_x_max) max(new_x_min, new_x_max)];
    new_YWorldLimits = [min(new_y_min, new_y_max) max(new_y_min, new_y_max)];

    verbose = 0;
     if verbose
        fprintf('%s%s%s%s%s\n', string(datetime(now,'ConvertFrom','datenum')), ' ', mfilename, ': ', 'Alignment adjustment:')
        fprintf('%0.1f%s%0.1f\n', x_celldim, ' ', y_celldim);
        fprintf('%i%s%i\n\n', x_offset_cells(1), ' ', y_offset_cells(1));
        %fprintf('%0.1f%s%0.1f\n\n', DEMs.refmat(3,1), ' ', DEMs.refmat(3,2));
        fprintf('%0.1f%s%0.1f\n', original_XWorldLimits(1), ' ', original_YWorldLimits(1));
        fprintf('%0.1f%s%0.1f\n\n\n', original_XWorldLimits(2), ' ', original_YWorldLimits(2));
        fprintf('%0.1f%s%0.1f\n', target_XWorldLimits(1), ' ', target_YWorldLimits(1));
        fprintf('%0.1f%s%0.1f\n\n\n', target_XWorldLimits(2), ' ', target_YWorldLimits(2));
        fprintf('%0.1f%s%0.1f\n', new_XWorldLimits(1), ' ', new_YWorldLimits(1));
        fprintf('%0.1f%s%0.1f\n\n\n', new_XWorldLimits(2), ' ', new_YWorldLimits(2));
    end

end
