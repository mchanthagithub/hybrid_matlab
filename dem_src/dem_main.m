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
% g = [0 0 0];

%% Material properties
rho = 1000;
% Spring constant
k = 5e6;
% Viscous dampener
gamma = 10.5;

n_dem = 200;
r_mean = 0.5;
r_sigma = 0.5*(0.15/3);
m_mean = (4/3)*pi*(r_mean^3)*rho;

t_col = pi*(2*k/m_mean - (gamma^2)/4)^(-1/2)
t_col_20 = t_col/20
rest_coeff = exp(-gamma*t_col/2)

%% Generate initial points
% Should have 4 state variables: q, r, v and m
% q is (n_DEM)x(3) in size
rngseed = 123;
[q,r] = generateInitialPointPositionsAndRadii(n_dem,r_mean,r_sigma,L,domain_min,domain_max,rngseed);
% Initialize random velocities
v = rand(n_dem,3);
% v = zeros(n_dem,3);
% v(:,3) = zeros(n_dem,1);
% q(1,:) = [0 0 5.501];
% q(2,:) = [1 0 4.501];
% q(1,:) = [0 0 5.5];
% q(2,:) = [2 0 5.5];
% r(1) = 0.5;
% r(2) = 0.5;
% v(1,:) = [1 0 0];
% v(2,:) = [-1 0 0];
mass = (4/3)*pi*(r.^3)*rho.*ones(1,n_dem);

%% Set any static planes (Only have axis-aligned planes)
static_planes(1).q = [0 0 0];
static_planes(1).n = [0 0 1];

static_planes(2).q = [10 0 0];
static_planes(2).n = [-1 0 0];

static_planes(3).q = [0 0 0];
static_planes(3).n = [1 0 0];

static_planes(4).q = [0 10 0];
static_planes(4).n = [0 -1 0];

static_planes(5).q = [0 0 0];
static_planes(5).n = [0 1 0];

static_planes(6).q = [0 0 20];
static_planes(6).n = [0 0 -1];

[dummy,n_static_planes] = size(static_planes);

q_save(1,:,:) = q;
v_save(1,:,:) = v;
ctr = 0;
save_per_sec = 25;
num_steps = floor(1.0/del_t);
save_freq = num_steps/save_per_sec;
save_ctr = 1;

%% Begin time-integration
for t = del_t:del_t:t_final
    F = zeros(n_dem,3);
    % Apply gravity
    for kk = 1:n_dem
        F(kk,:) = F(kk,:) + mass(kk)*g;
    end    
    % Look for sphere-sphere contacts
    contact_idx_sphere_sphere = zeros(0,0);
    contact_point_sphere_sphere = zeros(0,0);
    contact_normal_sphere_sphere = zeros(0,0);
    contact_pen_depth_sphere_sphere = zeros(0);
    contact_v_rel_sphere_sphere = zeros(0,0);
    contact_num_sphere_sphere = 0;
    for kk = 1:n_dem
        for jj = kk+1:n_dem
            diff = q(kk,:) - q(jj,:);
            dist = norm(diff);
            if(dist < r(jj) + r(kk))
                overlap = true;
                contact_num_sphere_sphere = contact_num_sphere_sphere + 1;
                contact_idx_sphere_sphere(contact_num_sphere_sphere,:) = [kk jj];
                contact_normal_sphere_sphere(contact_num_sphere_sphere,:) = diff./dist;
                contact_pen_depth_sphere_sphere(contact_num_sphere_sphere) = dist - (r(jj) + r(kk));
                contact_point_sphere_sphere(contact_num_sphere_sphere,:) = q(kk,:) + contact_normal_sphere_sphere(contact_num_sphere_sphere,:)*contact_pen_depth_sphere_sphere(contact_num_sphere_sphere)/2.0;
                contact_v_rel_sphere_sphere(contact_num_sphere_sphere,:) = v(kk,:) - v(jj,:);
            end
        end
    end
    
    % Look for sphere-plane contacts
    contact_idx_sphere_plane = zeros(0,0);
    contact_point_sphere_plane = zeros(0,0);
    contact_normal_sphere_plane = zeros(0,0);
    contact_pen_depth_sphere_plane = zeros(0);
    contact_v_rel_sphere_plane = zeros(0,0);
    contact_num_sphere_plane = 0;
    for kk = 1:n_dem
        for jj = 1:n_static_planes
            diff = q(kk,:) - static_planes(jj).q;
            dist = dot(static_planes(jj).n,diff);
            if(dist < r(kk))
                overlap = true;
                contact_num_sphere_plane = contact_num_sphere_plane + 1;
                contact_idx_sphere_plane(contact_num_sphere_plane,:) = [kk jj];
                contact_normal_sphere_plane(contact_num_sphere_plane,:) = static_planes(jj).n;
                contact_pen_depth_sphere_plane(contact_num_sphere_plane) = dist - r(kk);
%                 contact_pen_depth_sphere_plane(contact_num_sphere_plane)
                contact_point_sphere_plane(contact_num_sphere_plane,:) = q(kk,:) - contact_normal_sphere_plane(contact_num_sphere_plane,:)*contact_pen_depth_sphere_plane(contact_num_sphere_plane)/2.0;
                contact_v_rel_sphere_plane(contact_num_sphere_plane,:) = v(kk,:);
%                 contact_v_rel_sphere_plane(contact_num_sphere_plane,:)
            end
        end
    end
    
    
    % Loop over sphere-sphere contacts
    [num_contacts_sphere_sphere,] = size(contact_idx_sphere_sphere);
    for kk = 1:num_contacts_sphere_sphere
        v_normal = dot(contact_v_rel_sphere_sphere(kk,:),contact_normal_sphere_sphere(kk,:))*contact_normal_sphere_sphere(kk,:);
        indices = contact_idx_sphere_sphere(kk,:);
        idx1 = indices(1); idx2 = indices(2);
        m_eff = mass(idx1)*mass(idx2)/(mass(idx1)+mass(idx2));
        force_normal = -k * contact_pen_depth_sphere_sphere(kk) * contact_normal_sphere_sphere(kk,:) - ...
            gamma * m_eff * v_normal;
        F(idx1,:) = F(idx1,:) + force_normal;
        F(idx2,:) = F(idx2,:) - force_normal;
    end
    
    % Loop over sphere-plane contacts
    [num_contacts_sphere_plane,] = size(contact_idx_sphere_plane);
    for kk = 1:num_contacts_sphere_plane
%         kk
        v_normal = dot(contact_v_rel_sphere_plane(kk,:),contact_normal_sphere_plane(kk,:))*contact_normal_sphere_plane(kk,:);
%         -k * contact_pen_depth_sphere_plane(kk) * contact_normal_sphere_plane(kk,:)
%         -0.5 * gamma * v_normal
        force_normal = -k * contact_pen_depth_sphere_plane(kk) * contact_normal_sphere_plane(kk,:) - ...
            0.5 * gamma * mass(contact_idx_sphere_plane(kk,1)) * v_normal;
%         force_normal
        F(contact_idx_sphere_plane(kk),:) = F(contact_idx_sphere_plane(kk),:) + force_normal;
    end
%     F
    % Calculate new velocities (momentum update, forward euler)
    for kk = 1:n_dem
        v(kk,:) = v(kk,:) + del_t*(1/mass(kk))*F(kk,:);
    end    
%     v
    q = q + del_t*v;
%     q
%     if(ctr == 1)
%         return
%     end
    ctr = ctr+1;
    if(mod(ctr,save_freq) == 0)
        t
        save_ctr = save_ctr+1;
        q_save(save_ctr,:,:) = q;
        v_save(save_ctr,:,:) = v;
    end
end

for n = 1:save_ctr
   scatter3(q_save(n,:,1),q_save(n,:,2),q_save(n,:,3))
   zlim([-10,10])
   drawnow
end

for n = 1:save_ctr
    num_digits = numel(num2str(save_ctr));
    zero_string = '';
    for ii = 1:(num_digits - numel(num2str(n)))
        zero_string = strcat(zero_string,'0');
    end
    temp_string = strcat('dem_config_',zero_string);
    temp_string2 = strcat(temp_string,num2str(n));
    outputDEMVTK(strcat(temp_string2,'.vtu'),n_dem,squeeze(q_save(n,:,:)),squeeze(v_save(n,:,:)));
end

figure
hold on
plot(v_save(:,1,1))
plot(v_save(:,2,1))


