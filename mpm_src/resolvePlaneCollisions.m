function [ physics_grid ] = resolvePlaneCollisions(physics_grid,static_planes)
% Rasterize forces to grid
tol = 1e-14;
[dummy,num_static_planes] = size(static_planes);
total_node_num = prod(physics_grid.num_grid_nodes);
node_vels = physics_grid.rasterized_momentum_post_force./physics_grid.rasterized_mass;
% 8888888888888888
for node_num = 1:total_node_num
    if physics_grid.rasterized_mass(node_num) < tol
        continue
    end
    
    node_num_adjusted = node_num-1;
    node_idx = [mod(node_num_adjusted,physics_grid.num_grid_nodes(1)) ...
        mod(floor(node_num_adjusted/physics_grid.num_grid_nodes(1)), physics_grid.num_grid_nodes(2))...
        floor(node_num_adjusted/(physics_grid.num_grid_nodes(1)*physics_grid.num_grid_nodes(2)))];
%     node_idx
    
    node_coord = node_idx * physics_grid.delta + physics_grid.min;
    for plane_num = 1:num_static_planes
        dist = dot(static_planes(plane_num).n,node_coord -static_planes(plane_num).q);
        if(dist > 0.0)
            continue;
        end
        % sliding
        if(static_planes(plane_num).type == 0)
            vnormal = dot(static_planes(plane_num).n,node_vels(node_num,:));
            if(vnormal > 0.0)
                continue;
            end
            
            node_vels(node_num,:) = node_vels(node_num,:) - vnormal * static_planes(plane_num).n;
            physics_grid.rasterized_momentum_post_forces(node_num,:) = physics_grid.rasterized_mass(node_num) * node_vels(node_num,:);       
        % Sticking
        else
            physics_grid.rasterized_momentum_post_force(node_num,:) = 0;
        end        
    end
end
