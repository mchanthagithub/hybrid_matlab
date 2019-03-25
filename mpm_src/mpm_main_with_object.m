clear
close all
clc

L_grid = 5;
bodies(1).min = [2 2 2];
bodies(1).max = [3 3 3];
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
mpm_state = mpm_state(L_grid,bodies,static_planes,rho,K,nu,mu_s,mu_2,xi);

tic
%% Time integration
q_save(1,:,:) = mpm_state.mpm_points.q;
v_save(1,:,:) = mpm_state.mpm_points.vel;
ctr = 0;
save_per_sec = 25;
num_steps = floor(1.0/mpm_state.del_t);
save_freq = num_steps/save_per_sec;
save_ctr = 1;
tol = 1e-12;
num_saves_total = ceil(save_per_sec * mpm_state.t_final);

for t = mpm_state.del_t:mpm_state.del_t:mpm_state.t_final
    mpm_state = mpm_state.flowFirstPhase();
    mpm_state = mpm_state.flowSecondPhase();
    
    ctr = ctr+1;

    if(mod(ctr,save_freq) == 0)
        t
        save_ctr = save_ctr+1;
        q_save(save_ctr,:,:) = mpm_state.mpm_points.q;
        v_save(save_ctr,:,:) = mpm_state.mpm_points.vel;
        
        num_digits = numel(num2str(num_saves_total));
        zero_string = '';
        for ii = 1:(num_digits - numel(num2str(save_ctr)))
            zero_string = strcat(zero_string,'0');
        end
        temp_string = strcat('results/mpm_config_',zero_string);
        temp_string2 = strcat(temp_string,num2str(save_ctr));
        sigma_voigt(:,1) = mpm_state.mpm_points.sigma(:,1,1);
        sigma_voigt(:,2) = mpm_state.mpm_points.sigma(:,2,2);
        sigma_voigt(:,3) = mpm_state.mpm_points.sigma(:,3,3);
        sigma_voigt(:,4) = mpm_state.mpm_points.sigma(:,2,3);
        sigma_voigt(:,5) = mpm_state.mpm_points.sigma(:,3,1);
        sigma_voigt(:,6) = mpm_state.mpm_points.sigma(:,1,2);
        outputMPMVTK(strcat(temp_string2,'.vtu'),mpm_state.mpm_points.num_points,squeeze(q_save(save_ctr,:,:)),squeeze(v_save(save_ctr,:,:)),...
            sigma_voigt,mpm_state.mpm_points.gamma_bar_dot_p);
    end
end

toc