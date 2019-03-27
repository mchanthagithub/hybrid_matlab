%% Code to do 3d DEM
clc
clear
close all

%% Simulation parameters
g = [0 0 -9.81];

n_dem = 3000;
L = 10;
domain_min = [0,0,0]; %x,y,z
domain_max = [L,L,L]; %x,y,z
r_mean = 0.5;
r_sigma = 0.5*(0.15/3);
rho = 1000;
k_n = 5e6;
gamma_n = 10.5;
k_t = k_n/3.5;
gamma_t = gamma_n/2.0;
mu = 0.5;
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
[dummy,static_planes.num] = size(static_planes);
dem_state = dem_state(1,'results/old/dem_config_035.vtu','results/old/dem_contact_035.txt',n_dem,domain_min,domain_max,r_mean,r_sigma,...
    k_n,k_t,gamma_n,gamma_t,rho,mu,static_planes);
material_parameters = [dem_state.k_n dem_state.k_t dem_state.gamma_n dem_state.gamma_t dem_state.rho dem_state.mu];

q_save(1,:,:) = dem_state.q;
v_save(1,:,:) = dem_state.v;
ctr = 0;
save_per_sec = 50;
num_steps = floor(1.0/dem_state.del_t);
save_freq = num_steps/save_per_sec;
save_ctr = 1;
num_saves_total = ceil(save_per_sec * dem_state.t_final);

[state_file_string, contact_file_string] = makeDEMOutputFileStrings(num_saves_total,save_ctr,'results/dem_config_','results/dem_contact_');
outputDEMVTK(state_file_string,dem_state.n_dem,squeeze(q_save(save_ctr,:,:)),dem_state.r,squeeze(v_save(save_ctr,:,:)),material_parameters);
outputDEMContacts(contact_file_string,dem_state.contact_sphere_sphere_old,dem_state.contact_sphere_plane_old,dem_state.static_planes);

%% Begin time-integration
for t = dem_state.del_t:dem_state.del_t:dem_state.t_final
    dem_state = dem_state.flow();
    ctr = ctr+1;
    if(mod(ctr,save_freq) == 0)
        t
        dem_state.contact_sphere_sphere.num
        dem_state.contact_sphere_plane.num
        save_ctr = save_ctr+1;
        q_save(save_ctr,:,:) = dem_state.q;
        v_save(save_ctr,:,:) = dem_state.v;
        
        [state_file_string, contact_file_string] = makeDEMOutputFileStrings(num_saves_total,save_ctr,'results/dem_config_','results/dem_contact_');
        outputDEMVTK(state_file_string,dem_state.n_dem,squeeze(q_save(save_ctr,:,:)),dem_state.r,squeeze(v_save(save_ctr,:,:)),material_parameters);
        outputDEMContacts(contact_file_string,dem_state.contact_sphere_sphere_old,dem_state.contact_sphere_plane_old,dem_state.static_planes);
        fprintf('Done writing outputfiles\n');
    end
end



