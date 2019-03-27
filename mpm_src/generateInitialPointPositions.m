function [ q ] = generateInitialPointPositions( body, physics_grid )
%% Generate Points
[dummy,num_bodies] = size(body);

num_points = 0;
point_delta = physics_grid.delta / (physics_grid.points_per_dim);
for bdy_num = 1:num_bodies
    elem_min_idx = floor((body(bdy_num).min - physics_grid.min)./physics_grid.delta);
    elem_max_idx = ceil((body(bdy_num).max - physics_grid.min)./physics_grid.delta);
    elem_min_idx
    elem_max_idx
    
    num_points_prev = num_points;
    num_new_points = prod(elem_max_idx-elem_min_idx)*physics_grid.points_per_dim^3;
    num_points = num_points + num_new_points;
    
    q((num_points_prev+1):num_points,:) = zeros(num_new_points,3);
    
    coord_idx = num_points_prev+1;
    ctr = 0;
    for z_idx = elem_min_idx(3):(elem_max_idx(3)-1)
        for y_idx = elem_min_idx(2):(elem_max_idx(2)-1)
            for x_idx = elem_min_idx(1):(elem_max_idx(1)-1)
                base_coords = [x_idx y_idx z_idx]*physics_grid.delta + physics_grid.min + point_delta/2.0;
                
                for z_inner_idx = 1:physics_grid.points_per_dim
                    for y_inner_idx = 1:physics_grid.points_per_dim
                        for x_inner_idx = 1:physics_grid.points_per_dim
                            q(coord_idx,:) = base_coords + point_delta*[(x_inner_idx-1) (y_inner_idx-1) (z_inner_idx-1)];
                            coord_idx = coord_idx+1;
                        end
                    end
                end
            end
        end
    end
end

