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
k_n = 5e6;
% Viscous dampener
gamma_n = 10.5;

k_t = k_n/3.5;
gamma_t = gamma_n/2.0;
mu = 0.5;

n_dem = 3000;
r_mean = 0.5;
r_sigma = 0.5*(0.15/3);
m_mean = (4/3)*pi*(r_mean^3)*rho;

t_col = pi*(2*k_n/m_mean - (gamma_n^2)/4)^(-1/2)
t_col_20 = t_col/20
del_t
rest_coeff = exp(-gamma_n*t_col/2)

%% Generate initial points
% Should have 4 state variables: q, r, v and m
% q is (n_DEM)x(3) in size
rngseed = 123;
rng(rngseed);
[q,r] = generateInitialPointPositionsAndRadiiVoxel( n_dem,r_mean,r_sigma,domain_min,domain_max,rngseed );
[n_dem,dummy] = size(q);
n_dem
% Initialize random velocities
v = rand(n_dem,3);
[contact_sphere_sphere_old, contact_sphere_plane_old] = intializeContacts();

% [q,r,v,contact_sphere_sphere_old,contact_sphere_plane_old] = inputDEMStateFromFile('results/old/dem_config_070.vtu','results/old/dem_contacts_070.txt');
% [n_dem,dummy] = size(q);
% n_dem

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

[dummy,n_static_planes] = size(static_planes);

q_save(1,:,:) = q;
v_save(1,:,:) = v;
ctr = 0;
save_per_sec = 50;
num_steps = floor(1.0/del_t);
save_freq = num_steps/save_per_sec;
save_ctr = 1;
num_saves_total = ceil(save_per_sec * t_final);

[state_file_string, contact_file_string] = makeDEMOutputFileStrings(num_saves_total,save_ctr,'results/dem_config_','results/dem_contact_');
outputDEMVTK(state_file_string,n_dem,squeeze(q_save(save_ctr,:,:)),r,squeeze(v_save(save_ctr,:,:)));
outputDEMContacts(contact_file_string,contact_sphere_sphere_old,contact_sphere_plane_old);

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
    aabb_min = q - r';
    aabb_max = q + r';
    
    % Rasterize bounding boxes to grid
    grid = rasterizeAABBsToGrid(aabb_min,aabb_max,r_mean);
    
    % Find sphere_sphere contacts
    contact_sphere_sphere = findSphereSphereCollisionsUsingGrid( contact_sphere_sphere,contact_sphere_sphere_old, grid,q,r,v,del_t );
    
    % Find sphere_plane contacts
%     contact_sphere_plane = findSpherePlaneCollisionsBruteForce( contact_sphere_plane,contact_sphere_plane_old, static_planes,q,r,v,del_t );
    
    % Look for sphere-plane contacts
    for kk = 1:n_dem
        for jj = 1:n_static_planes
            diff = q(kk,:) - static_planes(jj).q;
            dist = dot(static_planes(jj).n,diff);
            if(dist < r(kk))
                contact_sphere_plane.num = contact_sphere_plane.num + 1;
                idx = [kk jj];
                normal = static_planes(jj).n;
                pen_depth = dist - r(kk);
                point = q(kk,:) - normal*pen_depth/2.0;
                v_rel = v(kk,:);
                
                contact_sphere_plane.idx(contact_sphere_plane.num,:) = idx;
                contact_sphere_plane.normal(contact_sphere_plane.num,:) = normal;
                contact_sphere_plane.pen_depth(contact_sphere_plane.num) = pen_depth;
%                 contact_sphere_plane_pen_depth_sphere_plane(contact_sphere_plane_num_sphere_plane)
                contact_sphere_plane.point(contact_sphere_plane.num,:) = point;
                contact_sphere_plane.v_rel(contact_sphere_plane.num,:) = v(kk,:);
%                 contact_v_rel_sphere_plane(contact_num_sphere_plane,:)
                contact_exists = 0;
                if(contact_sphere_plane_old.num > 0)
                    [contact_exists, index] = ismember([kk jj],contact_sphere_plane_old.idx,'rows');
                end
                if(contact_exists)
%                     index
                    delta_s = contact_sphere_plane_old.delta_s(index,:);
                    delta_s = delta_s + del_t*v_rel;
                    delta_s = delta_s - dot(normal,delta_s)*normal;
                    contact_sphere_plane.delta_s(contact_sphere_plane.num,:) = delta_s;
                else
                    delta_s = del_t * v_rel;
                    delta_s = delta_s - dot(normal,delta_s)*normal;
                    contact_sphere_plane.delta_s(contact_sphere_plane.num,:) = delta_s;
                end
            end
        end
    end
    
    
    % Loop over sphere-sphere contacts
    for kk = 1:contact_sphere_sphere.num
        indices = contact_sphere_sphere.idx(kk,:);
        idx1 = indices(1); idx2 = indices(2);
        m_eff = mass(idx1)*mass(idx2)/(mass(idx1)+mass(idx2));
        
        n = contact_sphere_sphere.normal(kk,:);
        delta_s = contact_sphere_sphere.delta_s(kk,:);
        pen_depth = contact_sphere_sphere.pen_depth(kk);
        vel = contact_sphere_sphere.v_rel(kk,:);
        vt = vel - dot(n,vel)*n;
        vn = vel - vt;
        
        normal_force = -k_n * pen_depth * n - gamma_n * m_eff * vn;
        friction_force = -k_t * delta_s - gamma_t * m_eff * vn;
        
        % Project into friction cone
        mu_fn = mu * norm(normal_force);
        ft = norm(friction_force);
        if(0.5 * m_eff * gamma_t * norm(vt) > mu_fn)
           contact_sphere_sphere.delta_s(kk,:) = 0;
           friction_force = -0.5 * m_eff * gamma_t * vt;
           friction_force = friction_force./norm(friction_force) * mu_fn;
        elseif(ft > mu_fn)
            a = k_t * k_t * norm(delta_s)^2;
            b = k_t * gamma_t * m_eff * dot(delta_s,vt);
            c = 0.25 * gamma_t * gamma_t * m_eff * m_eff * dot(vt,vt) - mu * mu * norm(normal_force)^2;
            dscr = b * b - 4.0 * a * c;
            dscr_sqrt = sqrt(dscr);
            root1 = secondRootOfQuadratic(a,b,c,dscr_sqrt);
            new_delta_s = root1 * delta_s;
            new_friction_force = -k_t * new_delta_s - 0.5 * m_eff * gamma_t * vt;
            friction_force = new_friction_force;
            contact_sphere_sphere.delta_s(kk,:) = new_delta_s;
        end
%         F(contact_sphere_sphere.idx(kk),:) = F(contact_sphere_sphere.idx(kk),:) + normal_force + friction_force;
        F(idx1,:) = F(idx1,:) + normal_force + friction_force;
        F(idx2,:) = F(idx2,:) - normal_force - friction_force;
    end
    
    % Loop over sphere-plane contacts
    for kk = 1:contact_sphere_plane.num
        m_eff = mass(contact_sphere_plane.idx(kk,1));
        n = contact_sphere_plane.normal(kk,:);
        delta_s = contact_sphere_plane.delta_s(kk,:);
        pen_depth = contact_sphere_plane.pen_depth(kk);
        vel = contact_sphere_plane.v_rel(kk,:);
        vt = vel - dot(n,vel)*n;
        vn = vel - vt;
        
        normal_force = -k_n * pen_depth * n - ...
            0.5 * gamma_n * m_eff * vn;
        
        friction_force = -k_t * delta_s - 0.5 * m_eff * gamma_t * vt;
        
        % Project into friction cone
        mu_fn = mu * norm(normal_force);
        ft = norm(friction_force);
        if(0.5 * m_eff * gamma_t * norm(vt) > mu_fn)
           contact_sphere_plane.delta_s(kk,:) = 0;
           friction_force = -0.5 * m_eff * gamma_t * vt;
           friction_force = friction_force./norm(friction_force) * mu_fn;
        elseif(ft > mu_fn)
            a = k_t * k_t * norm(delta_s)^2;
            b = k_t * gamma_t * m_eff * dot(delta_s,vt);
            c = 0.25 * gamma_t * gamma_t * m_eff * m_eff * dot(vt,vt) - mu * mu * norm(normal_force)^2;
            dscr = b * b - 4.0 * a * c;
            dscr_sqrt = sqrt(dscr);
            root1 = secondRootOfQuadratic(a,b,c,dscr_sqrt);
            new_delta_s = root1 * delta_s;
            new_friction_force = -k_t * new_delta_s - 0.5 * m_eff * gamma_t * vt;
            friction_force = new_friction_force;
            contact_sphere_plane.delta_s(kk,:) = new_delta_s;
        end
        F(contact_sphere_plane.idx(kk),:) = F(contact_sphere_plane.idx(kk),:) + normal_force + friction_force;
    end
    contact_sphere_sphere_old = contact_sphere_sphere;
    contact_sphere_plane_old = contact_sphere_plane;
    
    % Calculate new velocities (momentum update, forward euler)
    for kk = 1:n_dem
        v(kk,:) = v(kk,:) + del_t*(1/mass(kk))*F(kk,:);
    end    
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
        outputDEMVTK(state_file_string,n_dem,squeeze(q_save(save_ctr,:,:)),r,squeeze(v_save(save_ctr,:,:)));
        outputDEMContacts(contact_file_string,contact_sphere_sphere_old,contact_sphere_plane_old);
    end
end
toc
% for n = 1:save_ctr
%     num_digits = numel(num2str(save_ctr));
%     zero_string = '';
%     for ii = 1:(num_digits - numel(num2str(n)))
%         zero_string = strcat(zero_string,'0');
%     end
%     temp_string = strcat('dem_config_',zero_string);
%     temp_string2 = strcat(temp_string,num2str(n));
%     outputDEMVTK(strcat(temp_string2,'.vtu'),n_dem,squeeze(q_save(n,:,:)),squeeze(v_save(n,:,:)));
% end
% 
% figure
% hold on
% plot(v_save(:,1,1))
% plot(v_save(:,2,1))


