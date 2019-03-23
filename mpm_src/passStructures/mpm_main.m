%% Code to 3D MPM
clc
clear
close all

%% Simulation parameters
t_final = 0.4;
del_t = 0.0002;
g = [0 0 -9.81];
flip_weight = 1.0;

%% Grid parameters
L = 5;
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

%% Material properties
rho = 1000;
% % Bulk modulus
% K = 70e6;
% % Poisson ratio
% nu = 0.3;
% % Shear modulus
% G = 3*K*(1-2*nu)/(2*(1+nu));
points_per_dim = 2;

% E = 3*K*(1-2*nu);
material_properties.rho = 1000;
material_properties.K = 60e6;
material_properties.nu = 0.3;
material_properties.G = 3*material_properties.K*(1-2*material_properties.nu)/(2*(1+material_properties.nu));
material_properties.E = 3*material_properties.K*(1-2*material_properties.nu);
material_properties.rho_critical = material_properties.rho;

material_properties.mu_s = 0.3819;
material_properties.mu_2 = 0.6435;
material_properties.xi = 1.1233;


cfl_t = physics_grid.delta * sqrt(rho/material_properties.E)
del_t

%% Make bodies
body(1).min = [2 2 2];
body(1).max = [3 3 3];

% body(2).min = [2 5 1];
% body(2).max = [3 6 2];

%% Generate Points
[ q ] = generateInitialPointPositions( body, physics_grid );
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
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

%% Basis functions
basis_functions = createBasisFunctions;

%% Set any static planes (Only have axis-aligned planes)
% type = 0: sliding; type = 1:sticking
static_planes(1).q = [0 0 0];
static_planes(1).n = [0 0 1];
static_planes(1).type = 1;
% 
% static_planes(2).q = [10 0 0];
% static_planes(2).n = [-1 0 0];
% static_planes(2).type = 0;
% 
% static_planes(3).q = [0 0 0];
% static_planes(3).n = [1 0 0];
% static_planes(3).type = 0;
% 
% static_planes(4).q = [0 10 0];
% static_planes(4).n = [0 -1 0];
% static_planes(4).type = 0;
% 
% static_planes(5).q = [0 0 0];
% static_planes(5).n = [0 1 0];
% static_planes(5).type = 0;
% 
% static_planes(6).q = [0 0 20];
% static_planes(6).n = [0 0 -1];
% static_planes(6).type = 0;

[dummy,n_static_planes] = size(static_planes);

tic

%% Time integration
q_save(1,:,:) = mpm_points.q;
v_save(1,:,:) = mpm_points.vel;
ctr = 0;
save_per_sec = 25;
num_steps = floor(1.0/del_t);
save_freq = num_steps/save_per_sec;
save_ctr = 1;
tol = 1e-12;

num_saves_total = ceil(save_per_sec * t_final);

for t = del_t:del_t:t_final
    % Clear physics grid
    physics_grid.rasterized_mass(:) = 0;
    physics_grid.rasterized_momentum(:) = 0;
    physics_grid.rasterized_momentum_post_force(:) = 0;
    physics_grid.rasterized_forces(:) = 0;
    physics_grid.rasterized_velocity(:) = 0;
    physics_grid.rasterized_acceleration(:) = 0;    
    
    % Grid: Project mass to grid
    physics_grid = rasterizeMassToGrid(physics_grid,mpm_points,basis_functions);
    % Grid: Project momentum to grid
    physics_grid = rasterizeMomentumToGrid(physics_grid,mpm_points,basis_functions);
    
    % Points: Update point volume
    for pt_num = 1:mpm_points.num_points
        mpm_points.volume(pt_num) = mpm_points.volume(pt_num).*exp(del_t*trace(squeeze(mpm_points.vel_grad(pt_num,:,:))));
    end
%     physics_grid.rasterized_momentum(physics_grid.rasterized_momentum<0)
    % Points: Update stress
%     mpm_points = computeHypoelasticCauchyStressMuOfI(mpm_points,material_properties,del_t);
    mpm_points = computeHypoelasticCauchyStressLinearElastic(mpm_points,material_properties,del_t);
    % Grid: Project forces to grid
    physics_grid = computeForces(physics_grid,mpm_points,basis_functions,g);
    % Grid: Update grid momentum
    physics_grid.rasterized_momentum_post_force = physics_grid.rasterized_momentum + del_t * physics_grid.rasterized_forces;
%     physics_grid.rasterized_momentum_post_force(physics_grid.rasterized_momentum_post_force<0)
    % Grid: Resolve collisions with planes
    physics_grid = resolvePlaneCollisions(physics_grid,static_planes);
    
    % Grid: Update grid velocity
    physics_grid.rasterized_velocity = physics_grid.rasterized_momentum_post_force./physics_grid.rasterized_mass;
    physics_grid.rasterized_velocity(isnan(physics_grid.rasterized_velocity)) = 0;
    % Grid: Update grid acceleration
    physics_grid.rasterized_acceleration = (physics_grid.rasterized_momentum_post_force -physics_grid.rasterized_momentum)./(del_t*physics_grid.rasterized_mass);
    physics_grid.rasterized_acceleration(isnan(physics_grid.rasterized_acceleration)) = 0;
    
    % Points: Update velocity gradient
    mpm_points = computeVelocityGradient( physics_grid,mpm_points,basis_functions );
    
    % Points: Update velocities
    mpm_points = updatePointVelocities( physics_grid,mpm_points,basis_functions, del_t, flip_weight);
    
    % Points: Update positions
    mpm_points = updatePointPositions( physics_grid,mpm_points,basis_functions, del_t);
    
    % Points: Resolve plane collisions
    mpm_points = resolvePlanePointCollisions(mpm_points,static_planes);
    
    ctr = ctr+1;
%     if(ctr == 1)
%        mpm_points_new(ctr) = mpm_points;
%     elseif(ctr == 2)
%         mpm_points_new(ctr) = mpm_points;
%         return
%     end
%     
    if(mod(ctr,save_freq) == 0)
        t
        save_ctr = save_ctr+1;
        q_save(save_ctr,:,:) = mpm_points.q;
        v_save(save_ctr,:,:) = mpm_points.vel;
        
        num_digits = numel(num2str(num_saves_total));
        zero_string = '';
        for ii = 1:(num_digits - numel(num2str(save_ctr)))
            zero_string = strcat(zero_string,'0');
        end
        temp_string = strcat('results/mpm_config_',zero_string);
        temp_string2 = strcat(temp_string,num2str(save_ctr));
        sigma_voigt(:,1) = mpm_points.sigma(:,1,1);
        sigma_voigt(:,2) = mpm_points.sigma(:,2,2);
        sigma_voigt(:,3) = mpm_points.sigma(:,3,3);
        sigma_voigt(:,4) = mpm_points.sigma(:,2,3);
        sigma_voigt(:,5) = mpm_points.sigma(:,3,1);
        sigma_voigt(:,6) = mpm_points.sigma(:,1,2);
        outputMPMVTK(strcat(temp_string2,'.vtu'),mpm_points.num_points,squeeze(q_save(save_ctr,:,:)),squeeze(v_save(save_ctr,:,:)),...
            sigma_voigt,mpm_points.gamma_bar_dot_p);
    end
end


% for n = 1:save_ctr
%    scatter3(q_save(n,:,1),q_save(n,:,2),q_save(n,:,3))
%    zlim([0,10])
%    drawnow
% end

% for n = 1:save_ctr
%     num_digits = numel(num2str(save_ctr));
%     zero_string = '';
%     for ii = 1:(num_digits - numel(num2str(n)))
%         zero_string = strcat(zero_string,'0');
%     end
%     temp_string = strcat('mpm_config_',zero_string);
%     temp_string2 = strcat(temp_string,num2str(n));
%     outputMPMVTK(strcat(temp_string2,'.vtu'),mpm_points.num_points,squeeze(q_save(n,:,:)),squeeze(v_save(n,:,:)));
% end

% figure
% hold on
% plot(v_save(:,1,3))
toc