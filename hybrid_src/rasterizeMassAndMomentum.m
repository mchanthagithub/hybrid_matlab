function [rasterized_mass,rasterized_momentum,dem_coupled] = rasterizeMassAndMomentum( dem_state_q,dem_state_v,dem_state_mass,physics_grid_min, physics_grid_delta, physics_grid_num_grid_nodes, physics_grid_rasterized_mass,mpm_points_momentum,mpm_points_num_points,basis_functions )
    
    [n_dem,~] = size(dem_state_q);
    rasterized_mass = zeros(prod(physics_grid_num_grid_nodes),1);
    rasterized_momentum = zeros(prod(physics_grid_num_grid_nodes),3);
    dem_coupled = ones(n_dem,1);
    
    tol = 1e-14;
    for grain_num = 1:n_dem
        % Indices of lower left corner node of element the mpm point is in
        grid_idx = floor((dem_state_q(grain_num,:) - physics_grid_min)/physics_grid_delta);
        % Coordinates of lower left node
        min_node_coord = [grid_idx(1) grid_idx(2) grid_idx(3)]*physics_grid_delta + physics_grid_min;
        
        % Coordinates of q in a reference element (-1,-1,-1) to (1,1,1) assuming a cube element
        q_reference = 2*(dem_state_q(grain_num,:) - min_node_coord)/physics_grid_delta - 1;
        
        all_nodes_have_mpm_mass = 1;
        for z_idx = grid_idx(3):grid_idx(3)+1
            for y_idx = grid_idx(2):grid_idx(2)+1
                for x_idx = grid_idx(1):grid_idx(1)+1
                    grid_flat_idx = x_idx+1 + (y_idx)*physics_grid_num_grid_nodes(1) + (z_idx)*physics_grid_num_grid_nodes(1)*physics_grid_num_grid_nodes(2);
                    if(physics_grid_rasterized_mass(grid_flat_idx) < tol)
                        all_nodes_have_mpm_mass = 0;
                        break
                    end
                end
                if(all_nodes_have_mpm_mass == 0)
                    break
                end
            end
        end
        
        if(all_nodes_have_mpm_mass == 0)
            dem_coupled(grain_num) = 0;
            continue
        end
        
        for z_idx = grid_idx(3):grid_idx(3)+1
            for y_idx = grid_idx(2):grid_idx(2)+1
                for x_idx = grid_idx(1):grid_idx(1)+1
                    element_node_flat_idx = (x_idx-grid_idx(1)+1) + (y_idx-grid_idx(2))*2 + (z_idx-grid_idx(3))*4;
                                        
                    grid_flat_idx = x_idx+1 + (y_idx)*physics_grid_num_grid_nodes(1) + (z_idx)*physics_grid_num_grid_nodes(1)*physics_grid_num_grid_nodes(2);
                    
                    rasterized_mass(grid_flat_idx) = rasterized_mass(grid_flat_idx) + ...
                        basis_functions.N{element_node_flat_idx}(q_reference)*dem_state_mass(grain_num);
                    
                    rasterized_momentum(grid_flat_idx,:) = rasterized_momentum(grid_flat_idx,:) + ...
                        basis_functions.N{element_node_flat_idx}(q_reference)*dem_state_v(grain_num,:)*dem_state_mass(grain_num);
                end
            end
        end
        
        
    end
end

