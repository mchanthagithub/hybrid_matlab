classdef dem_state
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        q; % Position
        r; % radii
        v; % velocities
        contact_sphere_sphere;
        contact_sphere_plane;
        contact_sphere_sphere_old; % sphere-sphere contact cache
        contact_sphere_plane_old; % sphere-plane contact cache
        static_planes;
        n_dem;
        k_n; % normal spring stiffness coeff
        k_t; % tangential spring stifness coeff
        gamma_n; % normal damping coeff
        gamma_t; % tangential damping coeff
        rho; % density 
        mu; % friction coefficient
        F; % Forces
        mass;
        t_final = 10.0;
        del_t = 0.001;
        g = [0 0 -9.81];
    end
    
    methods
        function obj = dem_state(restart,state_filename,contact_filename,n_dem_in,domain_min,domain_max,r_mean,r_sigma,k_n,k_t,gamma_n,gamma_t,rho,mu,static_planes)
            if(restart)
                [q_read,r_read,v_read,material_parameters,contact_sphere_sphere_old_read, ...
                contact_sphere_plane_old_read,static_planes_read] = inputDEMStateFromFile(state_filename,contact_filename);
                [obj.n_dem,~] = size(q_read);
                obj.q = q_read;
                obj.r = r_read;
                obj.v = v_read;
                obj.k_n = material_parameters(1);
                obj.k_t = material_parameters(2);
                obj.gamma_n = material_parameters(3);
                obj.gamma_t = material_parameters(4);
                obj.rho = material_parameters(5);
                obj.mu = material_parameters(6);
                obj.contact_sphere_sphere_old = contact_sphere_sphere_old_read;
                obj.contact_sphere_plane_old = contact_sphere_plane_old_read;
                obj.static_planes = static_planes_read;
                obj.mass = (4/3)*pi*(obj.r.^3)*obj.rho.*ones(obj.n_dem,1);
                obj.F = zeros(obj.n_dem,3);
            else
                rngseed = 123;
                rng(rngseed);
                obj.k_n = k_n; obj.k_t = k_t; obj.gamma_n = gamma_n; obj.gamma_t = gamma_t; 
                obj.rho = rho; obj.mu = mu;
                [obj.q,obj.r] = generateInitialPointPositionsAndRadiiVoxel( n_dem_in,r_mean,r_sigma,domain_min,domain_max,rngseed );
                [obj.n_dem,~] = size(obj.q);
                % Initialize random velocities
                obj.v = rand(obj.n_dem,3);
                obj.F = zeros(obj.n_dem,3);
                [obj.contact_sphere_sphere_old, obj.contact_sphere_plane_old] = intializeContacts();
                [obj.contact_sphere_sphere, obj.contact_sphere_plane] = intializeContacts();
                obj.static_planes = static_planes;
                obj.mass = (4/3)*pi*(obj.r.^3)*obj.rho.*ones(obj.n_dem,1);
            end
            
        end
        
%         function obj = dem_state_new(n_dem,domain_min,domain_max,r_mean,r_sigma,k_n,k_t,gamma_n,gamma_t,rho,mu,static_planes)
%             rngseed = 123;
%             rng(rngseed);
%             obj.k_n = k_n; obj.k_t = k_t; obj.gamma_n = gamma_n; obj.gamma_t = gamma_t; 
%             obj.rho = rho; obj.mu = mu;
%             [obj.q,obj.r] = generateInitialPointPositionsAndRadiiVoxel( n_dem,r_mean,r_sigma,domain_min,domain_max,rngseed );
%             [obj.n_dem,~] = size(obj.q);
%             % Initialize random velocities
%             obj.v = rand(n_dem,3);
%             obj.F = zeros(n_dem,3);
%             [obj.contact_sphere_sphere_old, obj.contact_sphere_plane_old] = intializeContacts();
%             [obj.contact_sphere_sphere, obj.contact_sphere_plane] = intializeContacts();
%             obj.static_planes = static_planes;
%             obj.mass = (4/3)*pi*(obj.r.^3)*obj.rho.*ones(obj.n_dem,1);
%         end
%         
        function obj = flow(obj)
            r_mean = mean(obj.r);
            obj.F = zeros(obj.n_dem,3);
            % Apply gravity
            for kk = 1:obj.n_dem
                obj.F(kk,:) = obj.F(kk,:) + obj.mass(kk)*obj.g;
            end

            [obj.contact_sphere_sphere, obj.contact_sphere_plane] = intializeContacts();

            % Set up bounding boxes
            aabb_min = obj.q - obj.r;
            aabb_max = obj.q + obj.r;

            % Rasterize bounding boxes to grid
            grid = rasterizeAABBsToGrid(aabb_min,aabb_max,r_mean);

            % Find sphere_sphere contacts
            obj.contact_sphere_sphere = findSphereSphereCollisionsUsingGrid( obj.contact_sphere_sphere,obj.contact_sphere_sphere_old, grid,obj.q,obj.r,obj.v,obj.del_t );

            % Find sphere_plane contacts
            obj.contact_sphere_plane = findSpherePlaneCollisionsBruteForce( obj.contact_sphere_plane,obj.contact_sphere_plane_old, obj.static_planes,obj.q,obj.r,obj.v,obj.del_t );

            % Calculate forces from sphere-sphere contacts
            [obj.F,obj.contact_sphere_sphere] = calculateSphereSphereForces(obj.F,obj.contact_sphere_sphere,obj.mass,obj.k_n,obj.k_t,obj.gamma_n,obj.gamma_t,obj.mu);

            % Calculate forces from sphere-plane contacts
            [obj.F,obj.contact_sphere_plane] = calculateSpherePlaneForces(obj.F,obj.contact_sphere_plane,obj.mass,obj.k_n,obj.k_t,obj.gamma_n,obj.gamma_t,obj.mu);

            % Cache old contacts
            obj.contact_sphere_sphere_old = obj.contact_sphere_sphere;
            obj.contact_sphere_plane_old = obj.contact_sphere_plane;

            % Calculate new velocities (momentum update, forward euler)
            for kk = 1:obj.n_dem
                obj.v(kk,:) = obj.v(kk,:) + obj.del_t*(1/obj.mass(kk))*obj.F(kk,:);
            end

            % Calculate new positions
            obj.q = obj.q + obj.del_t*obj.v;
        end
        
        % Grabs dem grains based on idx
        function obj_out = dem_state_splice(obj,dem_forward_index_map,dem_reverse_index_map)
            obj_out = obj;
            obj_out.q = obj.q(dem_forward_index_map,:);
            [obj_out.n_dem,~] = size(obj.q(dem_forward_index_map,:));
            obj_out.r = obj.r(dem_forward_index_map);
            obj_out.v = obj.v(dem_forward_index_map,:);
            obj_out.F = obj.F(dem_forward_index_map,:);
            obj_out.mass = obj.mass(dem_forward_index_map);
            
            new_contact_ctr = 0;
            % Adjust indices of cached sphere sphere collisions
            for ii = 1:obj.contact_sphere_sphere_old.num
                idx1 = obj.contact_sphere_sphere_old.idx(ii,1);
                idx2 = obj.contact_sphere_sphere_old.idx(ii,2);
                new_idx1 = dem_reverse_index_map(idx1);
                new_idx2 = dem_reverse_index_map(idx2);
                if(new_idx1 == 0 || new_idx2 == 0)
                   continue; 
                end
                new_contact_ctr = new_contact_ctr + 1;
                obj_out.contact_sphere_sphere_old.idx(new_contact_ctr,:) = [new_idx1 new_idx2];
            end
            obj_out.contact_sphere_sphere_old.idx = obj_out.contact_sphere_sphere_old.idx(1:new_contact_ctr,:);
            obj_out.contact_sphere_sphere_old.num = new_contact_ctr;
            
            new_contact_ctr = 0;
            % Adjust indices of cached sphere plane collisions
            for ii = 1:obj.contact_sphere_plane_old.num
                idx1 = obj.contact_sphere_plane_old.idx(ii,1);
                idx2 = obj.contact_sphere_plane_old.idx(ii,2);
                new_idx1 = dem_reverse_index_map(idx1);
                %idx2 is for the plane so don't adjust that
                if(new_idx1 == 0)
                   continue; 
                end
                new_contact_ctr = new_contact_ctr + 1;
                obj_out.contact_sphere_plane_old.idx(new_contact_ctr,:) = [new_idx1 idx2];
            end
            obj_out.contact_sphere_plane_old.idx = obj_out.contact_sphere_plane_old.idx(1:new_contact_ctr,:);
            obj_out.contact_sphere_plane_old.num = new_contact_ctr;
        end
        
    end
end

