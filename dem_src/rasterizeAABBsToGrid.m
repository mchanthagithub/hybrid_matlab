function [ grid ] = rasterizeAABBsToGrid( aabb_min, aabb_max, r_mean)
[num_grains,dummy] = size(aabb_min);
delta_h = r_mean * 3.0;
grid_min = [realmax realmax realmax];
grid_max = [realmin realmin realmin];
for grain_num = 1:num_grains
    grid_min = min(grid_min,aabb_min(grain_num,:));
    grid_max = max(grid_max,aabb_max(grain_num,:));
end


grid.num_grid_elements = ceil((grid_max - grid_min)/delta_h);
grid.total_num_grid_elements = prod(grid.num_grid_elements);
grid.collisions = cell(grid.total_num_grid_elements,1);
aabb_element_idx_min = floor((aabb_min - grid_min)./delta_h);
aabb_element_idx_min(aabb_element_idx_min < 0) = 0;
aabb_element_idx_max = ceil((aabb_max - grid_min)./delta_h);

for grain_num = 1:num_grains
   for kk = aabb_element_idx_min(grain_num,3):(aabb_element_idx_max(grain_num,3)-1)
       for jj = aabb_element_idx_min(grain_num,2):(aabb_element_idx_max(grain_num,2)-1)
           for ii = aabb_element_idx_min(grain_num,1):(aabb_element_idx_max(grain_num,1)-1)
               global_element_flat_idx = ii+1 + jj*grid.num_grid_elements(1) + kk*grid.num_grid_elements(1)*grid.num_grid_elements(2);
               grid.collisions{global_element_flat_idx} = [grid.collisions{global_element_flat_idx} grain_num];
           end
       end
   end
end

