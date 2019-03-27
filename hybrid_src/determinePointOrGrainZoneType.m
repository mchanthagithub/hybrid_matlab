function [pt_zone_type] = determinePointOrGrainZoneType(q,physics_grid_min, physics_grid_delta,...
    physics_grid_num_grid_elements,element_zone_type)
    [num_total_pts,~] = size(q);
    pt_zone_type = zeros(num_total_pts,1);
    
    total_grid_elements = prod(physics_grid_num_grid_elements);
    for pt_num = 1:num_total_pts
        % Indices of lower left corner node of element the mpm point is in
        grid_idx = floor((q(pt_num,:) - physics_grid_min)/physics_grid_delta);
        global_element_num = grid_idx(1)+1 + grid_idx(2)*physics_grid_num_grid_elements(1) + ...
            grid_idx(3)*physics_grid_num_grid_elements(1)*physics_grid_num_grid_elements(2);
        
        pt_zone_type(pt_num) = element_zone_type(global_element_num);
    end
    

end

