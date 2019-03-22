function [ mpm_points ] = computeHypoelasticCauchyStressLinearElastic( mpm_points,material_properties,dt )
function [dev] = deviator(A)
    dev = A - (trace(A)/3)*eye(3);
end

E = material_properties.E;
nu = material_properties.nu;
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
    
    mpm_points.sigma(pnt_idx,:,:) = sigma_trial;
end
end

