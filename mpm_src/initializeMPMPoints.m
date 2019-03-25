function [mpm_points] = initializeMPMPoints(bodies,physics_grid,material_properties)
[ q ] = generateInitialPointPositions( bodies, physics_grid );
mpm_points.q = q;
[num_points,dummy] = size(q);
mpm_points.num_points = num_points;
mpm_points.momentum = zeros(size(q));
mpm_points.vel = zeros(size(q));
mpm_points.rho = material_properties.rho*ones(num_points,1);
mpm_points.volume = (physics_grid.delta^3)/(physics_grid.points_per_dim^3) * ones(num_points,1);
mpm_points.mass = mpm_points.rho .* mpm_points.volume;
mpm_points.sigma = zeros(num_points,3,3);
mpm_points.vel_grad = zeros(num_points,3,3);
mpm_points.gamma_bar_dot_p = zeros(num_points,1);
scatter3(q(:,1),q(:,2),q(:,3))
end

