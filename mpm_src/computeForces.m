function [ physics_grid ] = computeForces( physics_grid,mpm_points,basis_functions,g )
    % Rasterize forces to grid
    for pt_num = 1:mpm_points.num_points
        
        % Indices of lower left corner node of element the mpm point is in
        grid_idx = floor((mpm_points.q(pt_num,:) - physics_grid.min)/physics_grid.delta);
        % Coordinates of lower left node
        min_node_coord = [grid_idx(1) grid_idx(2) grid_idx(3)]*physics_grid.delta + physics_grid.min;
        
        % Coordinates of q in a reference element (-1,-1,-1) to (1,1,1) assuming a cube element
        q_reference = 2*(mpm_points.q(pt_num,:) - min_node_coord)/physics_grid.delta - 1;
        
        basis_grad = zeros(3,1);
        sigma = squeeze(mpm_points.sigma(pt_num,:,:));
        for z_idx = grid_idx(3):grid_idx(3)+1
            for y_idx = grid_idx(2):grid_idx(2)+1
                for x_idx = grid_idx(1):grid_idx(1)+1
                    element_node_flat_idx = (x_idx-grid_idx(1)+1) + (y_idx-grid_idx(2))*2 + (z_idx-grid_idx(3))*4;
                                        
                    grid_flat_idx = x_idx + (y_idx-1)*physics_grid.num_grid_nodes(1) + (z_idx-1)*physics_grid.num_grid_nodes(1)*physics_grid.num_grid_nodes(2);
                    
                    basis_grad(1) = basis_functions.dN{element_node_flat_idx,1}(q_reference);
                    basis_grad(2) = basis_functions.dN{element_node_flat_idx,2}(q_reference);
                    basis_grad(3) = basis_functions.dN{element_node_flat_idx,3}(q_reference);
                    
                    internal_force = -mpm_points.volume(pt_num) * (sigma*basis_grad);
                    
                    physics_grid.rasterized_forces(grid_flat_idx,:) = physics_grid.rasterized_forces(grid_flat_idx,:) + internal_force';
                end
            end
        end
    end
    physics_grid.rasterized_forces = physics_grid.rasterized_forces + physics_grid.rasterized_mass*g;
end

