%% Code to do 3d DEM
clc
clear
close all

%% Simulation parameters
t_final = 10.0;
del_t = 0.001;
L = 10;
domain_min = [0,0,0]; %x,y,z
domain_max = [L,L,L]; %x,y,z
g = [0 0 -9.81];

%% Generate initial points
% Should have 4 state variables: q, r, v and m
% q is (n_DEM)x(3) in size
% rngseed = 123;
% rng(rngseed);
% n_dem = 3000;
% r_mean = 0.5;
% r_sigma = 0.5*(0.15/3);
% [q,r] = generateInitialPointPositionsAndRadiiVoxel( n_dem,r_mean,r_sigma,domain_min,domain_max,rngseed );
% [n_dem,dummy] = size(q);
% n_dem
% % Initialize random velocities
% v = rand(n_dem,3);
% [contact_sphere_sphere_old, contact_sphere_plane_old] = intializeContacts();
% % Material properties
% rho = 1000;
% k_n = 5e6;
% gamma_n = 10.5;
% k_t = k_n/3.5;
% gamma_t = gamma_n/2.0;
% mu = 0.5;
% material_parameters = [k_n k_t gamma_n gamma_t rho mu];
% 
% m_mean = (4/3)*pi*(r_mean^3)*rho;
%% Set any static planes (Only have axis-aligned planes)
% static_planes(1).q = [0 0 0];
% static_planes(1).n = [0 0 1];
% 
% static_planes(2).q = [10 0 0];
% static_planes(2).n = [-1 0 0];
% 
% static_planes(3).q = [0 0 0];
% static_planes(3).n = [1 0 0];
% 
% static_planes(4).q = [0 10 0];
% static_planes(4).n = [0 -1 0];
% 
% static_planes(5).q = [0 0 0];
% static_planes(5).n = [0 1 0];
% 
% [dummy,static_planes.num] = size(static_planes);


%% Read in points
[q,r,v,material_parameters,contact_sphere_sphere_old,contact_sphere_plane_old,static_planes] = inputDEMStateFromFile('results/old/dem_config_040.vtu','results/old/dem_contact_040.txt');
[n_dem,dummy] = size(q);
n_dem
k_n = material_parameters(1)
k_t = material_parameters(2)
gamma_n = material_parameters(3)
gamma_t = material_parameters(4)
rho = material_parameters(5)
mu = material_parameters(6)

%%
mass = (4/3)*pi*(r.^3)*rho.*ones(n_dem,1);
m_mean = mean(mass);
t_col = pi*(2*k_n/m_mean - (gamma_n^2)/4)^(-1/2)
t_col_20 = t_col/20
del_t
rest_coeff = exp(-gamma_n*t_col/2)


q_save(1,:,:) = q;
v_save(1,:,:) = v;
ctr = 0;
save_per_sec = 50;
num_steps = floor(1.0/del_t);
save_freq = num_steps/save_per_sec;
save_ctr = 1;
num_saves_total = ceil(save_per_sec * t_final);

[state_file_string, contact_file_string] = makeDEMOutputFileStrings(num_saves_total,save_ctr,'results/dem_config_','results/dem_contact_');
outputDEMVTK(state_file_string,n_dem,squeeze(q_save(save_ctr,:,:)),r,squeeze(v_save(save_ctr,:,:)),material_parameters);
outputDEMContacts(contact_file_string,contact_sphere_sphere_old,contact_sphere_plane_old,static_planes);

tic
r_mean = mean(r);
%% Begin time-integration
for t = del_t:del_t:t_final
    F = zeros(n_dem,3);
    % Apply gravity
    for kk = 1:n_dem
        F(kk,:) = F(kk,:) + mass(kk)*g;
    end
    
    [contact_sphere_sphere, contact_sphere_plane] = intializeContacts();
    
    % Set up bounding boxes
    aabb_min = q - r;
    aabb_max = q + r;
    
    % Rasterize bounding boxes to grid
    grid = rasterizeAABBsToGrid(aabb_min,aabb_max,r_mean);
    
    % Find sphere_sphere contacts
    contact_sphere_sphere = findSphereSphereCollisionsUsingGrid( contact_sphere_sphere,contact_sphere_sphere_old, grid,q,r,v,del_t );
    
    % Find sphere_plane contacts
    contact_sphere_plane = findSpherePlaneCollisionsBruteForce( contact_sphere_plane,contact_sphere_plane_old, static_planes,q,r,v,del_t );
    
    % Calculate forces from sphere-sphere contacts
    [F,contact_sphere_sphere] = calculateSphereSphereForces(F,contact_sphere_sphere,mass,k_n,k_t,gamma_n,gamma_t,mu);
    
    % Calculate forces from sphere-plane contacts
    [F,contact_sphere_plane] = calculateSpherePlaneForces(F,contact_sphere_plane,mass,k_n,k_t,gamma_n,gamma_t,mu);
    
    % Cache old contacts
    contact_sphere_sphere_old = contact_sphere_sphere;
    contact_sphere_plane_old = contact_sphere_plane;
    
    % Calculate new velocities (momentum update, forward euler)
    for kk = 1:n_dem
        v(kk,:) = v(kk,:) + del_t*(1/mass(kk))*F(kk,:);
    end
    
    % Calculate new positions
    q = q + del_t*v;

    ctr = ctr+1;
    if(mod(ctr,save_freq) == 0)
        t
        contact_sphere_sphere.num
        contact_sphere_plane.num
        save_ctr = save_ctr+1;
        q_save(save_ctr,:,:) = q;
        v_save(save_ctr,:,:) = v;
        
        [state_file_string, contact_file_string] = makeDEMOutputFileStrings(num_saves_total,save_ctr,'results/dem_config_','results/dem_contact_');
        outputDEMVTK(state_file_string,n_dem,squeeze(q_save(save_ctr,:,:)),r,squeeze(v_save(save_ctr,:,:)),material_parameters);
        outputDEMContacts(contact_file_string,contact_sphere_sphere_old,contact_sphere_plane_old,static_planes);
    end
end
toc

