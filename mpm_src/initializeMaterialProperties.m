function [material_properties] = initializeMaterialProperties(rho,K,nu,mu_s,mu_2,xi)
material_properties.rho = 1000;
material_properties.K = 60e6;
material_properties.nu = 0.3;
material_properties.G = 3*material_properties.K*(1-2*material_properties.nu)/(2*(1+material_properties.nu));
material_properties.E = 3*material_properties.K*(1-2*material_properties.nu);
material_properties.rho_critical = material_properties.rho;

material_properties.mu_s = 0.3819;
material_properties.mu_2 = 0.6435;
material_properties.xi = 1.1233;
end

