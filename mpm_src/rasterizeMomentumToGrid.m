function [ physics_grid ] = rasterizeMomentumToGrid( physics_grid,mpm_points,basis_functions )
    % Rasterize momentum to grid
    for pt_num = 1:mpm_points.num_points
        
        % Indices of lower left corner node of element the mpm point is in
        grid_idx = floor((mpm_points.q(pt_num,:) - physics_grid.min)/physics_grid.delta);
        % Coordinates of lower left node
        min_node_coord = [grid_idx(1) grid_idx(2) grid_idx(3)]*physics_grid.delta + physics_grid.min;
        
        % Coordinates of q in a reference element (-1,-1,-1) to (1,1,1) assuming a cube element
        q_reference = 2*(mpm_points.q(pt_num,:) - min_node_coord)/physics_grid.delta - 1;
        
        for z_idx = grid_idx(3):grid_idx(3)+1
            for y_idx = grid_idx(2):grid_idx(2)+1
                for x_idx = grid_idx(1):grid_idx(1)+1
                    element_node_flat_idx = (x_idx-grid_idx(1)+1) + (y_idx-grid_idx(2))*2 + (z_idx-grid_idx(3))*4;
                                        
                    grid_flat_idx = x_idx + (y_idx-1)*physics_grid.num_grid_nodes(1) + (z_idx-1)*physics_grid.num_grid_nodes(1)*physics_grid.num_grid_nodes(2);
                    physics_grid.rasterized_momentum(grid_flat_idx,:) = physics_grid.rasterized_momentum(grid_flat_idx,:) + ...
                        basis_functions.N{element_node_flat_idx}(q_reference)*mpm_points.momentum(pt_num,:);
                end
            end
        end
    end
    
end

