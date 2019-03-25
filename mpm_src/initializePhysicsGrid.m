function [physics_grid] = initializePhysicsGrid(L)
physics_grid.min = [0,0,0];
physics_grid.max = [L,L,L];
physics_grid.delta = L/10;
physics_grid.num_grid_nodes = ceil((physics_grid.max - physics_grid.min)./physics_grid.delta);
physics_grid.num_grid_elements = physics_grid.num_grid_nodes-1;
physics_grid.points_per_dim = 2;
physics_grid.rasterized_mass = zeros(prod(physics_grid.num_grid_nodes),1);
physics_grid.rasterized_momentum = zeros(prod(physics_grid.num_grid_nodes),3);
physics_grid.rasterized_momentum_post_force = zeros(prod(physics_grid.num_grid_nodes),3);
physics_grid.rasterized_forces = zeros(prod(physics_grid.num_grid_nodes),3);
physics_grid.rasterized_velocity = zeros(prod(physics_grid.num_grid_nodes),3);
physics_grid.rasterized_acceleration = zeros(prod(physics_grid.num_grid_nodes),3);
end
