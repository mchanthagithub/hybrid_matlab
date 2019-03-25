function [mpm_state] = flowFirstPhase(mpm_state)
    % Clear physics grid
    mpm_state.physics_grid = clearPhysicsGrid(mpm_state.physics_grid);
    
    % Grid: Project mass to grid
    mpm_state.physics_grid.rasterized_mass = rasterizeMassToGrid( mpm_state.physics_grid.min, mpm_state.physics_grid.delta, ...
        mpm_state.physics_grid.num_grid_nodes, mpm_state.mpm_points.q,mpm_state.mpm_points.mass,mpm_state.mpm_points.num_points,basis_functions );
    
    % Grid: Project momentum to grid
    mpm_state.physics_grid.rasterized_momentum = rasterizeMomentumToGrid( mpm_state.physics_grid.min, mpm_state.physics_grid.delta, ...
        mpm_state.physics_grid.num_grid_nodes, mpm_state.mpm_points.q,mpm_state.mpm_points.momentum,mpm_state.mpm_points.num_points,basis_functions );
    
    % Points: Update point volume
    for pt_num = 1:mpm_state.mpm_points.num_points
        mpm_state.mpm_points.volume(pt_num) = mpm_state.mpm_points.volume(pt_num).*exp(del_t*trace(squeeze(mpm_state.mpm_points.vel_grad(pt_num,:,:))));
    end
    
    % Points: Update stress
%     mpm_state.mpm_points = computeHypoelasticCauchyStressMuOfI(mpm_state.mpm_points,material_properties,del_t);
    mpm_state.mpm_points = computeHypoelasticCauchyStressLinearElastic(mpm_state.mpm_points,material_properties,del_t);
    
    % Grid: Project forces to grid
    mpm_state.physics_grid.rasterized_forces = computeForces( mpm_state.physics_grid.min, mpm_state.physics_grid.delta, mpm_state.physics_grid.num_grid_nodes, ...
        mpm_state.physics_grid.rasterized_mass, mpm_state.mpm_points.q,mpm_state.mpm_points.sigma,mpm_state.mpm_points.volume,mpm_state.mpm_points.num_points,basis_functions,g );
    
    % Grid: Update grid momentum
    mpm_state.physics_grid.rasterized_momentum_post_force = mpm_state.physics_grid.rasterized_momentum + del_t * mpm_state.physics_grid.rasterized_forces;
%     physics_grid.rasterized_momentum_post_force(physics_grid.rasterized_momentum_post_force<0)
    
% Grid: Resolve collisions with planes
    mpm_state.physics_grid = resolvePlaneCollisions(mpm_state.physics_grid,static_planes);
end

