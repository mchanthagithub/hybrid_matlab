function [ vel_new, momentum_new ] = updatePointVelocities( physics_grid_min, physics_grid_delta,...
    physics_grid_num_grid_nodes, physics_grid_rasterized_velocity,physics_grid_rasterized_acceleration, ...
    mpm_points_q,mpm_points_vel,mpm_points_mass,mpm_points_num_points,basis_functions, del_t, flip_weight)
    
    vel_new = zeros(mpm_points_num_points,3);
    momentum_new = zeros(mpm_points_num_points,3);
% Rasterize momentum to grid
    for pt_num = 1:mpm_points_num_points
        
        % Indices of lower left corner node of element the mpm point is in
        grid_idx = floor((mpm_points_q(pt_num,:) - physics_grid_min)/physics_grid_delta);
        % Coordinates of lower left node
        min_node_coord = [grid_idx(1) grid_idx(2) grid_idx(3)]*physics_grid_delta + physics_grid_min;
        
        % Coordinates of q in a reference element (-1,-1,-1) to (1,1,1) assuming a cube element
        q_reference = 2*(mpm_points_q(pt_num,:) - min_node_coord)/physics_grid_delta - 1;
        
        v_pic = zeros(1,3);
        accel = zeros(1,3);
        for z_idx = grid_idx(3):grid_idx(3)+1
            for y_idx = grid_idx(2):grid_idx(2)+1
                for x_idx = grid_idx(1):grid_idx(1)+1
                    element_node_flat_idx = (x_idx-grid_idx(1)+1) + (y_idx-grid_idx(2))*2 + (z_idx-grid_idx(3))*4;

                    grid_flat_idx = x_idx+1 + (y_idx)*physics_grid_num_grid_nodes(1) + (z_idx)*physics_grid_num_grid_nodes(1)*physics_grid_num_grid_nodes(2);
                    
                    basis_weight = basis_functions.N{element_node_flat_idx}(q_reference);

                    v_pic = v_pic + basis_weight*physics_grid_rasterized_velocity(grid_flat_idx,:);
                    accel = accel + basis_weight*physics_grid_rasterized_acceleration(grid_flat_idx,:);
                end
            end
        end
        v_flip = mpm_points_vel(pt_num,:) + accel*del_t;
        vel_new(pt_num,:) = (1.0 - flip_weight) * v_pic + flip_weight*v_flip;
        momentum_new(pt_num,:) = mpm_points_vel(pt_num,:) * mpm_points_mass(pt_num);
    end
    
end

