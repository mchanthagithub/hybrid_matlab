function [ rasterized_mass ] = rasterizeMassToGrid( physics_grid_min, physics_grid_delta, physics_grid_num_grid_nodes, mpm_points_q,mpm_points_mass,mpm_points_num_points,basis_functions )
    rasterized_mass = zeros(prod(physics_grid_num_grid_nodes),1);
% Rasterize mass to grid
    for pt_num = 1:mpm_points_num_points
        
        % Indices of lower left corner node of element the mpm point is in
        grid_idx = floor((mpm_points_q(pt_num,:) - physics_grid_min)/physics_grid_delta);
%         grid_idx(grid_idx<0) = 0;
        % Coordinates of lower left node
        min_node_coord = [grid_idx(1) grid_idx(2) grid_idx(3)]*physics_grid_delta + physics_grid_min;
        
        % Coordinates of q in a reference element (-1,-1,-1) to (1,1,1) assuming a cube element
        q_reference = 2*(mpm_points_q(pt_num,:) - min_node_coord)/physics_grid_delta - 1;
        
        for z_idx = grid_idx(3):grid_idx(3)+1
            for y_idx = grid_idx(2):grid_idx(2)+1
                for x_idx = grid_idx(1):grid_idx(1)+1
                    element_node_flat_idx = (x_idx-grid_idx(1)+1) + (y_idx-grid_idx(2))*2 + (z_idx-grid_idx(3))*4;
                                        
                    grid_flat_idx = x_idx+1 + (y_idx)*physics_grid_num_grid_nodes(1) + (z_idx)*physics_grid_num_grid_nodes(1)*physics_grid_num_grid_nodes(2);
                    if(grid_flat_idx < 1 || grid_flat_idx > prod(physics_grid_num_grid_nodes))
                        grid_idx
                    end
                    rasterized_mass(grid_flat_idx) = rasterized_mass(grid_flat_idx) + ...
                        basis_functions.N{element_node_flat_idx}(q_reference)*mpm_points_mass(pt_num);
                end
            end
        end
    end
    
end

