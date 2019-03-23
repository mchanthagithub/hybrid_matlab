function [connectivity, global_coordinates] = createElements(Lx,Ly,Lz,del_x,del_y,del_z)
if(mod(Lx,del_x) ~= 0 || mod(Ly,del_y) ~= 0 || mod(Lz,del_z) ~= 0)
    fprintf('Not an integer partition of element lengths')
    return
end

num_elem_x = Lx/del_x;
num_elem_y = Ly/del_y;
num_elem_z = Lz/del_z;
num_node_x = num_elem_x+1;
num_node_y = num_elem_y+1;
num_node_z = num_elem_z+1;

num_elements = num_elem_x*num_elem_y*num_elem_z;
num_grid_nodes = (num_elem_x+1)*(num_elem_y+1)*(num_elem_z+1);

base_coords = [0 0 0];
global_coordinates = zeros(num_elements,3);
connectivity = zeros(num_elements,8);

for kk = 1:num_elem_z
    for jj = 1:num_elem_y
        for ii = 1:num_elem_x
            global_elem_flat_idx = ii + (jj-1)*num_elem_x + (kk-1)*num_elem_x*num_elem_y;
            
            global_node_flat_idx_bottom_left = ii + (jj-1)*num_node_x + (kk-1)*num_node_x*num_node_y;
            
            connectivity(global_elem_flat_idx,1) = global_node_flat_idx_bottom_left;
            connectivity(global_elem_flat_idx,2) = global_node_flat_idx_bottom_left+1;
            connectivity(global_elem_flat_idx,3) = global_node_flat_idx_bottom_left + num_node_x+1;
            connectivity(global_elem_flat_idx,4) = global_node_flat_idx_bottom_left + num_node_x;
            
            connectivity(global_elem_flat_idx,5:8) = connectivity(global_elem_flat_idx,1:4) + num_node_x*num_node_z;
        end
    end
end

