function [ vel_grad ] = computeVelocityGradient( physics_grid_min, physics_grid_delta,...
    physics_grid_num_grid_nodes, physics_grid_rasterized_velocity,mpm_points_q,mpm_points_num_points,basis_functions )
  vel_grad = zeros(mpm_points_num_points,3,3);  
% Rasterize momentum to grid
    for pt_num = 1:mpm_points_num_points
        
        % Indices of lower left corner node of element the mpm point is in
        grid_idx = floor((mpm_points_q(pt_num,:) - physics_grid_min)/physics_grid_delta);
        % Coordinates of lower left node
        min_node_coord = [grid_idx(1) grid_idx(2) grid_idx(3)]*physics_grid_delta + physics_grid_min;
        
        % Coordinates of q in a reference element (-1,-1,-1) to (1,1,1) assuming a cube element
        q_reference = 2*(mpm_points_q(pt_num,:) - min_node_coord)/physics_grid_delta - 1;
        
        basis_grad = zeros(3,1);
        vel_grad_temp = zeros(3,3);
        for z_idx = grid_idx(3):grid_idx(3)+1
            for y_idx = grid_idx(2):grid_idx(2)+1
                for x_idx = grid_idx(1):grid_idx(1)+1
                    element_node_flat_idx = (x_idx-grid_idx(1)+1) + (y_idx-grid_idx(2))*2 + (z_idx-grid_idx(3))*4;
                    
                    grid_flat_idx = x_idx+1 + (y_idx)*physics_grid_num_grid_nodes(1) + (z_idx)*physics_grid_num_grid_nodes(1)*physics_grid_num_grid_nodes(2);
                    
                    basis_grad(1) = basis_functions.dN{element_node_flat_idx,1}(q_reference);
                    basis_grad(2) = basis_functions.dN{element_node_flat_idx,2}(q_reference);
                    basis_grad(3) = basis_functions.dN{element_node_flat_idx,3}(q_reference);
                    
                    vel_grad_temp = vel_grad_temp + basis_grad*physics_grid_rasterized_velocity(grid_flat_idx,:);
                end
            end
        end
        vel_grad(pt_num,:,:) = vel_grad_temp;
%         if(sum(vel_grad(vel_grad ~= 0))~=0)
%             vel_grad(vel_grad ~= 0)
%         end
    end
    
end

