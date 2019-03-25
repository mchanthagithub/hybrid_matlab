classdef mpm_state
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mpm_points;
        physics_grid;
        material_properties;
        static_planes;
        basis_functions;
        t_final = 0.4;
        del_t = 0.0002;
        g = [0 0 -9.81];
        flip_weight = 1.0;
    end
    
    methods
        function obj = mpm_state(L_grid,bodies,static_planes,rho,K,nu,mu_s,mu_2,xi)
            obj.material_properties = initializeMaterialProperties(rho,K,nu,mu_s,mu_2,xi);
            obj.physics_grid = initializePhysicsGrid(L_grid);
            obj.mpm_points = initializeMPMPoints(bodies,obj.physics_grid,obj.material_properties);
            obj.static_planes = static_planes;
            obj.basis_functions = createBasisFunctions();
        end
        
        function obj = flowFirstPhase(obj)
            % Clear physics grid
            obj.physics_grid = clearPhysicsGrid(obj.physics_grid);

            % Grid: Project mass to grid
            obj.physics_grid.rasterized_mass = rasterizeMassToGrid( obj.physics_grid.min, obj.physics_grid.delta, ...
                obj.physics_grid.num_grid_nodes, obj.mpm_points.q,obj.mpm_points.mass,obj.mpm_points.num_points,obj.basis_functions );

            % Grid: Project momentum to grid
            obj.physics_grid.rasterized_momentum = rasterizeMomentumToGrid( obj.physics_grid.min, obj.physics_grid.delta, ...
                obj.physics_grid.num_grid_nodes, obj.mpm_points.q,obj.mpm_points.momentum,obj.mpm_points.num_points,obj.basis_functions );

            % Points: Update point volume
            for pt_num = 1:obj.mpm_points.num_points
                obj.mpm_points.volume(pt_num) = obj.mpm_points.volume(pt_num).*exp(obj.del_t*trace(squeeze(obj.mpm_points.vel_grad(pt_num,:,:))));
            end

            % Points: Update stress
        %     obj.mpm_points = computeHypoelasticCauchyStressMuOfI(obj.mpm_points,material_properties,del_t);
            obj.mpm_points = computeHypoelasticCauchyStressLinearElastic(obj.mpm_points,obj.material_properties,obj.del_t);

            % Grid: Project forces to grid
            obj.physics_grid.rasterized_forces = computeForces( obj.physics_grid.min, obj.physics_grid.delta, obj.physics_grid.num_grid_nodes, ...
                obj.physics_grid.rasterized_mass, obj.mpm_points.q,obj.mpm_points.sigma,obj.mpm_points.volume,obj.mpm_points.num_points,obj.basis_functions,obj.g );

            % Grid: Update grid momentum
            obj.physics_grid.rasterized_momentum_post_force = obj.physics_grid.rasterized_momentum + obj.del_t * obj.physics_grid.rasterized_forces;
        %     physics_grid.rasterized_momentum_post_force(physics_grid.rasterized_momentum_post_force<0)

            % Grid: Resolve collisions with planes
            obj.physics_grid = resolvePlaneCollisions(obj.physics_grid,obj.static_planes);
        end
        
        function obj = flowSecondPhase(obj)
            % Grid: Update grid velocity
            obj.physics_grid.rasterized_velocity = obj.physics_grid.rasterized_momentum_post_force./obj.physics_grid.rasterized_mass;
            obj.physics_grid.rasterized_velocity(isnan(obj.physics_grid.rasterized_velocity)) = 0;
            % Grid: Update grid acceleration
            obj.physics_grid.rasterized_acceleration = (obj.physics_grid.rasterized_momentum_post_force -obj.physics_grid.rasterized_momentum)./(obj.del_t*obj.physics_grid.rasterized_mass);
            obj.physics_grid.rasterized_acceleration(isnan(obj.physics_grid.rasterized_acceleration)) = 0;

            % Points: Update velocity gradient
            obj.mpm_points.vel_grad = computeVelocityGradient( obj.physics_grid.min, obj.physics_grid.delta,...
                obj.physics_grid.num_grid_nodes, obj.physics_grid.rasterized_velocity,obj.mpm_points.q,obj.mpm_points.num_points,obj.basis_functions );

            % Points: Update velocities
            [obj.mpm_points.vel,obj.mpm_points.momentum] = updatePointVelocities( obj.physics_grid.min, obj.physics_grid.delta,...
                obj.physics_grid.num_grid_nodes, obj.physics_grid.rasterized_velocity,obj.physics_grid.rasterized_acceleration, ...
                obj.mpm_points.q,obj.mpm_points.vel,obj.mpm_points.mass,obj.mpm_points.num_points,obj.basis_functions, obj.del_t, obj.flip_weight);

            % Points: Update positions
            obj.mpm_points.q = updatePointPositions( obj.physics_grid.min, obj.physics_grid.delta,...
                obj.physics_grid.num_grid_nodes, obj.physics_grid.rasterized_velocity, ...
                obj.mpm_points.q,obj.mpm_points.num_points,obj.basis_functions, obj.del_t);

            % Points: Resolve plane collisions
            obj.mpm_points = resolvePlanePointCollisions(obj.mpm_points,obj.static_planes);
        end
    end
end

