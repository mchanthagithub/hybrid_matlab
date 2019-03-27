function [mpm_new_rasterized_momentum,dem_v_after_coupling,lambda,map,logical_map,dem_coupled] = hybridCoupling(mpm_state,dem_state)
    hybrid_grid = initializePhysicsGrid(mpm_state.physics_grid.min,mpm_state.physics_grid.max,mpm_state.physics_grid.delta);
    
%     hybrid_grid.rasterized_mass = rasterizeMassToGrid( mpm_state.physics_grid.min, ...
%         mpm_state.physics_grid.delta, mpm_state.physics_grid.num_grid_nodes,...
%         dem_state.q,dem_state.mass,dem_state.n_dem,mpm_state.basis_functions );
%     
%     dem_momentum = dem_state.v.*dem_state.mass;
%     hybrid_grid.rasterized_momentum = rasterizeMomentumToGrid( mpm_state.physics_grid.min, ...
%         mpm_state.physics_grid.delta, mpm_state.physics_grid.num_grid_nodes,...
%         dem_state.q,dem_momentum,dem_state.n_dem,mpm_state.basis_functions );
    
    [hybrid_grid.rasterized_mass, hybrid_grid.rasterized_momentum,dem_coupled] = rasterizeMassAndMomentum(...
        dem_state.q, dem_state.v, dem_state.mass, mpm_state.physics_grid.min, mpm_state.physics_grid.delta, ...
        mpm_state.physics_grid.num_grid_nodes, mpm_state.physics_grid.rasterized_mass,...
        mpm_state.mpm_points.momentum,mpm_state.mpm_points.num_points,mpm_state.basis_functions);

    m_mpm = mpm_state.physics_grid.rasterized_mass;
    m_dem = hybrid_grid.rasterized_mass;
    mass_sum = m_mpm + m_dem;
%     mass_sum = hybrid_grid.rasterized_mass + mpm_state.physics_grid.rasterized_mass;
    map = find(m_mpm > 0);
    logical_map = m_mpm > 0;
    mass_sum = mass_sum(logical_map);
    m_mpm = m_mpm(logical_map);
    m_dem = m_dem(logical_map);
    
    
    lambda = (1./mass_sum) .* (mpm_state.physics_grid.rasterized_momentum_post_force(logical_map,:).*m_dem - ...
        hybrid_grid.rasterized_momentum(logical_map,:).*m_mpm);
    
    mpm_new_rasterized_momentum =  mpm_state.physics_grid.rasterized_momentum_post_force;
    mpm_new_rasterized_momentum(map,:) = mpm_state.physics_grid.rasterized_momentum_post_force(map,:) - lambda;
%     mpm_state.physics_grid.rasterized_momentum_post_force(map,:) = mpm_state.physics_grid.rasterized_momentum_post_force(map,:) + lambda;
    
    dem_grid_momentum_after = hybrid_grid.rasterized_momentum(map,:) + lambda;
    
    dem_grid_v_after = zeros(prod(hybrid_grid.num_grid_nodes),3);
    dem_grid_a_after = zeros(prod(hybrid_grid.num_grid_nodes),3);
    
    dem_grid_v_after(map,:) = dem_grid_momentum_after ./ m_dem;
    dem_grid_a_after(map,:) = lambda ./ m_dem;
    
    dem_v_before_coupling = dem_state.v;
    
    [dem_v_after_coupling, ~ ] = updatePointVelocitiesHybrid( mpm_state.physics_grid.min, mpm_state.physics_grid.delta,...
        mpm_state.physics_grid.num_grid_nodes, dem_grid_v_after,dem_grid_a_after, ...
        dem_state.q,dem_state.v,dem_state.mass,dem_state.n_dem,mpm_state.basis_functions, dem_state.del_t, mpm_state.flip_weight,...
        dem_coupled);
    
    dem_v_after_coupling(~dem_coupled,:) = dem_v_before_coupling(~dem_coupled,:);
    
%     updatePointVelocities( physics_grid_min, physics_grid_delta,...
%     physics_grid_num_grid_nodes, physics_grid_rasterized_velocity,physics_grid_rasterized_acceleration, ...
%     mpm_points_q,mpm_points_vel,mpm_points_mass,mpm_points_num_points,basis_functions, del_t, flip_weight,dem_coupled)
end

