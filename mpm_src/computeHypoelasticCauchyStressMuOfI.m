function [ mpm_points ] = computeHypoelasticCauchyStressMuOfI( mpm_points,material_properties,dt )
function [dev] = deviator(A)
    dev = A - (trace(A)/3)*eye(3);
end

E = material_properties.E;
nu = material_properties.nu;
% Following algorithm in Dunatunga et al (2015)
iden = eye(3);
for pnt_idx = 1:mpm_points.num_points
    %std::cout<<"pntIdx: "<<pnt_idx<<" ========================================="<<std::endl;
    mpm_points.rho(pnt_idx) = mpm_points.mass(pnt_idx)/mpm_points.volume(pnt_idx);
    % Symmetric part of velocity gradient; stretch tensor
    D = 0.5*(squeeze(mpm_points.vel_grad(pnt_idx,:,:))+squeeze(mpm_points.vel_grad(pnt_idx,:,:))');
    % Skew part of velocity gradient, spin tensor
    W = 0.5*(squeeze(mpm_points.vel_grad(pnt_idx,:,:))-squeeze(mpm_points.vel_grad(pnt_idx,:,:))');
    % Objective correction for Jaummann rate
    G = -squeeze(mpm_points.sigma(pnt_idx,:,:))*W + W*squeeze(mpm_points.sigma(pnt_idx,:,:));
    % Trial stress tensor
    sigma_trial = squeeze(mpm_points.sigma(pnt_idx,:,:)) + dt*( (E/(1.0+nu))*(D + (nu/(1.0-2.0*nu))*trace(D)*iden ) + G);
    % Trial equivalent shear stress
    tau_bar_trial =  sqrt(0.5 * sum(sum(deviator(sigma_trial).*deviator(sigma_trial))) );
    % Trial pressure (note same as end pressure due to plastic incompressibility
    p_trial = -(1.0/3.0)*trace(sigma_trial);
    if(mpm_points.rho(pnt_idx) < material_properties.rho_critical)
      % In extension
      mpm_points.sigma(pnt_idx,:,:) = 0;
      mpm_points.gamma_bar_dot_p(pnt_idx) = (1.0/(material_properties.G*dt))*(tau_bar_trial);
    else
      if(p_trial <= 0.0)
        % In extension
        mpm_points.sigma(pnt_idx,:,:) = 0;
        mpm_points.gamma_bar_dot_p(pnt_idx) = (1.0/(material_properties.G*dt))*(tau_bar_trial);
      else
        S0 = material_properties.mu_s*p_trial;
        if(tau_bar_trial <= S0)
          % Elastic
          mpm_points.sigma(pnt_idx,:,:) = sigma_trial;
          mpm_points.gamma_bar_dot_p(pnt_idx) = 0.0;
        else
          % Plastic flow
          S2 = material_properties.mu_2*p_trial;
          alpha = material_properties.xi*material_properties.shear_modulus*dt*sqrt(p_trial);
          B = S2 + tau_bar_trial + alpha;
          H = S2*tau_bar_trial + S0*alpha;
          tau_bar_new = 2.0*H/(B + sqrt(B*B - 4.0*H));
          mpm_points.gamma_bar_dot_p(pnt_idx) = (1.0/(material_properties.G*dt))*(tau_bar_trial-tau_bar_new);
          mpm_points.sigma(pnt_idx,:,:) = (tau_bar_new/tau_bar_trial)*deviator(sigma_trial) - p_trial*iden;
        end
      end
    end
end
end

