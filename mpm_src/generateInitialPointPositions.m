function [ q ] = generateInitialPointPositions( num_bodies,body )
num_points = 0;
point_delta = grid_delta / (points_per_dim+1);
for bdy_num = 1:num_bodies
    elem_min_idx = floor((body(bdy_num).min - grid_min)./grid_delta);
    elem_max_idx = ceil((body(bdy_num).max - grid_min)./grid_delta);
   
    num_points_prev = num_points;
    num_new_points = prod(elem_max_idx+1-elem_min_idx)*points_per_dim^3;
    num_points = num_points + num_new_points;
    
    q((num_points_prev+1):num_points,:) = zeros(num_new_points,3);
    
    coord_idx = num_points_prev+1;
    ctr = 0;
    for z_idx = elem_min_idx(3):elem_max_idx(3)
        for y_idx = elem_min_idx(2):elem_max_idx(2)
            for x_idx = elem_min_idx(1):elem_max_idx(1)
                base_coords = [x_idx-1 y_idx-1 z_idx-1]*grid_delta + grid_min;
                
                for z_inner_idx = 1:points_per_dim
                    for y_inner_idx = 1:points_per_dim
                        for x_inner_idx = 1:points_per_dim
                            q(coord_idx,:) = base_coords + point_delta*[x_inner_idx y_inner_idx z_inner_idx];
                            coord_idx = coord_idx+1;
                        end
                    end
                end
            end
        end
    end
end

end

