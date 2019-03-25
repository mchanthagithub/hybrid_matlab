function [physics_grid] = clearPhysicsGrid(physics_grid)
    physics_grid.rasterized_mass(:) = 0;
    physics_grid.rasterized_momentum(:) = 0;
    physics_grid.rasterized_momentum_post_force(:) = 0;
    physics_grid.rasterized_forces(:) = 0;
    physics_grid.rasterized_velocity(:) = 0;
    physics_grid.rasterized_acceleration(:) = 0;   
end

