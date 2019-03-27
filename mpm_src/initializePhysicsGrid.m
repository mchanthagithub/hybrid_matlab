function [physics_grid] = initializePhysicsGrid(grid_min,grid_max,delta_h)
physics_grid.min = grid_min;
physics_grid.max = grid_max;
physics_grid.delta = delta_h;
physics_grid.num_grid_nodes = floor((physics_grid.max - physics_grid.min)./physics_grid.delta)+1;
physics_grid.num_grid_elements = physics_grid.num_grid_nodes-1;
physics_grid.points_per_dim = 2;
physics_grid.rasterized_mass = zeros(prod(physics_grid.num_grid_nodes),1);
physics_grid.rasterized_momentum = zeros(prod(physics_grid.num_grid_nodes),3);
physics_grid.rasterized_momentum_post_force = zeros(prod(physics_grid.num_grid_nodes),3);
physics_grid.rasterized_forces = zeros(prod(physics_grid.num_grid_nodes),3);
physics_grid.rasterized_velocity = zeros(prod(physics_grid.num_grid_nodes),3);
physics_grid.rasterized_acceleration = zeros(prod(physics_grid.num_grid_nodes),3);
end
