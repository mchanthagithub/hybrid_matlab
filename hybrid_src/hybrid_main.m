clear
close all
clc

addpath('../mpm_src')
addpath('../dem_src')

% Set up DEM
dem_state = dem_state(1,'results/old/dem_config_458.vtu','results/old/dem_contact_458.txt');
dem_material_parameters = [dem_state.k_n dem_state.k_t dem_state.gamma_n dem_state.gamma_t dem_state.rho dem_state.mu];

% Set up MPM
bodies(1).min = [0 0 0];
bodies(1).max = [6 6 84];
points_per_dim = 2;
rho = 1000;
K = 60e6;
nu = 0.3;
mu_s = 0.3819;
mu_2 = 0.6435;
xi = 1.1233;
static_planes(1).q = [0 0 0];
static_planes(1).n = [0 0 1];
static_planes(1).type = 1;

static_planes(2).q = [6 0 0];
static_planes(2).n = [-1 0 0];
static_planes(2).type = 0;

static_planes(3).q = [0 0 0];
static_planes(3).n = [1 0 0];
static_planes(3).type = 0;

static_planes(4).q = [0 6 0];
static_planes(4).n = [0 -1 0];
static_planes(4).type = 0;

static_planes(5).q = [0 0 0];
static_planes(5).n = [0 1 0];
static_planes(5).type = 0;

grid_min = [0 0 0];
grid_max = [6 6 90];
delta_h = 6;
mpm_state = mpm_state(grid_min,grid_max,delta_h,bodies,static_planes,rho,K,nu,mu_s,mu_2,xi);

del_t = min(mpm_state.del_t,dem_state.del_t);
t_final = max(mpm_state.t_final,dem_state.t_final);
mpm_state.del_t = del_t; dem_state.del_t = del_t;
mpm_state.t_final = t_final; dem_state.t_final = t_final;


ctr = 0;
save_per_sec = 200;
num_steps = floor(1.0/dem_state.del_t);
save_freq = num_steps/save_per_sec;
save_ctr = 1;
num_saves_total = ceil(save_per_sec * dem_state.t_final);

del_t
t_final
dem_time = 0;
mpm_time = 0;
hybrid_coupling_time = 0;
enrichment_time = 0;
enrichment_ctr = 100;
enrichment_freq = 100;

% Set zone types; 0 == continuum, 1 == discrete, 2 == hybrid
zone_type = ones(prod(mpm_state.physics_grid.num_grid_elements),1);
zone_type(3:4) = 2;
zone_type(5:6) = 0;
zone_type(7:8) = 2;
dem_grain_zone_type_old = zeros(dem_state.n_dem,1);
mpm_pt_zone_type_old = zeros(mpm_state.mpm_points.num_points,1);

dem_grain_zone_type_old = determinePointOrGrainZoneType(dem_state.q,mpm_state.physics_grid.min, mpm_state.physics_grid.delta,...
            mpm_state.physics_grid.num_grid_elements,zone_type);
mpm_pt_zone_type_old = determinePointOrGrainZoneType(mpm_state.mpm_points.q,mpm_state.physics_grid.min, mpm_state.physics_grid.delta,...
            mpm_state.physics_grid.num_grid_elements,zone_type);

old_dem_index_list = 1:dem_state.n_dem;
dem_reverse_index_map = old_dem_index_list'.*(dem_grain_zone_type_old~=0);
dem_forward_index_map = find(dem_grain_zone_type_old~=0);
[num_new_dem,~] = size(dem_forward_index_map);
dem_reverse_index_map(dem_reverse_index_map > 0) = 1:num_new_dem;
dem_state = dem_state_splice(dem_state,dem_forward_index_map,dem_reverse_index_map);

old_mpm_index_list = 1:mpm_state.mpm_points.num_points;
mpm_reverse_index_map = old_mpm_index_list'.*(mpm_pt_zone_type_old~=1);
mpm_forward_index_map = find(mpm_pt_zone_type_old~=1);
[num_new_mpm,~] = size(mpm_forward_index_map);
mpm_reverse_index_map(mpm_reverse_index_map > 0) = 1:num_new_mpm;
mpm_state = mpm_state_splice(mpm_state,mpm_forward_index_map);

q_save(1,:,:) = dem_state.q;
v_save(1,:,:) = dem_state.v;
[state_file_string, contact_file_string] = makeDEMOutputFileStrings(num_saves_total,save_ctr,'results/dem_config_','results/dem_contact_');
outputDEMVTK(state_file_string,dem_state.n_dem,dem_state.q,dem_state.r,dem_state.v,dem_material_parameters);
outputDEMContacts(contact_file_string,dem_state.contact_sphere_sphere_old,dem_state.contact_sphere_plane_old,dem_state.static_planes);

[mpm_file_string, ~] = makeDEMOutputFileStrings(num_saves_total,save_ctr,'results/mpm_config_','results/mpm_contact_');
outputMPMVTK(mpm_file_string,mpm_state.mpm_points.num_points,mpm_state.mpm_points.q,mpm_state.mpm_points.vel,...
    mpm_state.mpm_points.sigma,mpm_state.mpm_points.gamma_bar_dot_p);
        

%% Begin time-integration
for t = dem_state.del_t:dem_state.del_t:dem_state.t_final
    
    tic
    % Enrichment/determine zone types
    if(mod(enrichment_ctr,enrichment_freq)==0)
        enrichment_ctr = 1;
        dem_grain_zone_type = zeros(dem_state.n_dem,1);
        mpm_pt_zone_type = zeros(mpm_state.mpm_points.num_points,1);

        dem_grain_zone_type = determinePointOrGrainZoneType(dem_state.q,mpm_state.physics_grid.min, mpm_state.physics_grid.delta,...
                    mpm_state.physics_grid.num_grid_elements,zone_type);
        mpm_pt_zone_type = determinePointOrGrainZoneType(mpm_state.mpm_points.q,mpm_state.physics_grid.min, mpm_state.physics_grid.delta,...
                    mpm_state.physics_grid.num_grid_elements,zone_type);
        
        old_dem_index_list = 1:dem_state.n_dem;
        dem_reverse_index_map = old_dem_index_list'.*(dem_grain_zone_type~=0);
        dem_forward_index_map = find(dem_grain_zone_type~=0);
        [num_new_dem,~] = size(dem_forward_index_map);
        dem_reverse_index_map(dem_reverse_index_map > 0) = 1:num_new_dem;
        dem_state = dem_state_splice(dem_state,dem_forward_index_map,dem_reverse_index_map);

        old_mpm_index_list = 1:mpm_state.mpm_points.num_points;
        mpm_reverse_index_map = old_mpm_index_list'.*(mpm_pt_zone_type~=1);
        mpm_forward_index_map = find(mpm_pt_zone_type~=1);
        [num_new_mpm,~] = size(mpm_forward_index_map);
        mpm_reverse_index_map(mpm_reverse_index_map > 0) = 1:num_new_mpm;
        mpm_state = mpm_state_splice(mpm_state,mpm_forward_index_map);
    end
    enrichment_time = enrichment_time + toc;
    
    % DEM flow
    tic
    dem_state = dem_state.flow();
    dem_time = dem_time + toc;
    
    % MPM first flow
    tic
    mpm_state = mpm_state.flowFirstPhase();
    mpm_time = mpm_time + toc;

    % Hybrid coupling
    tic
    dem_v_before_coupling = dem_state.v;
    mpm_mom_before_coupling = mpm_state.physics_grid.rasterized_momentum_post_force;
    [mpm_new_rasterized_momentum,dem_v_after_coupling,lambda,map,logical_map,dem_coupled] = hybridCoupling(mpm_state,dem_state);
    mpm_state.physics_grid.rasterized_momentum_post_force = mpm_new_rasterized_momentum;
    dem_state.v = dem_v_after_coupling;
    diff_dem = dem_state.v - dem_v_before_coupling;
    diff_mpm = mpm_state.physics_grid.rasterized_momentum_post_force - mpm_mom_before_coupling;
    dem_state.q = dem_state.q + dem_state.del_t * (dem_state.v - dem_v_before_coupling);
    hybrid_coupling_time = hybrid_coupling_time + toc;
    
%     return
    % MPM second flow
    tic
    mpm_state = mpm_state.flowSecondPhase();
    mpm_time = mpm_time + toc;
    
    
    ctr = ctr+1;
    enrichment_ctr = enrichment_ctr+1;
    if(mod(ctr,save_freq) == 0)
        fprintf('Current t: %g dem time: %g mpm time: %g hybrid coupling time: %g enrichment time: %g\n',...
            t,dem_time,mpm_time,hybrid_coupling_time,enrichment_time);
%         fprintf('S-S contacts: %d S-P contacts: %d\n',dem_state.contact_sphere_sphere.num,dem_state.contact_sphere_plane.num);
        dem_time = 0;
        mpm_time = 0;
        hybrid_coupling_time = 0;
        enrichment_time = 0;
        save_ctr = save_ctr+1;
%         q_save(save_ctr,:,:) = dem_state.q;
%         v_save(save_ctr,:,:) = dem_state.v;
        
        % Output DEM
        [state_file_string, contact_file_string] = makeDEMOutputFileStrings(num_saves_total,save_ctr,'results/dem_config_','results/dem_contact_');
        outputDEMVTK(state_file_string,dem_state.n_dem,dem_state.q,dem_state.r,dem_state.v,dem_material_parameters);
        outputDEMContacts(contact_file_string,dem_state.contact_sphere_sphere_old,dem_state.contact_sphere_plane_old,dem_state.static_planes);
        
        % Output MPM
        [mpm_file_string, ~] = makeDEMOutputFileStrings(num_saves_total,save_ctr,'results/mpm_config_','results/mpm_contact_');
        outputMPMVTK(mpm_file_string,mpm_state.mpm_points.num_points,mpm_state.mpm_points.q,mpm_state.mpm_points.vel,...
            mpm_state.mpm_points.sigma,mpm_state.mpm_points.gamma_bar_dot_p);
        
        fprintf('Done writing outputfiles\n\n');
    end
end
