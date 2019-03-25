function [q,r] = generateInitialPointPositionsAndRadiiVoxel( n_dem,r_mean,r_sigma,domain_min,domain_max,rngseed )
rng(rngseed)
q = zeros(0,0);
r = normrnd(r_mean,r_sigma,1,n_dem);
r_max = max(r);
d_max = r_max*2.2;
domain_diff = domain_max - domain_min;
num_grains_in_dir = floor(domain_diff./(d_max));
ctr = 1;

z_coord = r_max + domain_min(3);
while(ctr <= n_dem)
    y_coord = r_max + domain_min(2);
    for jj = 1:num_grains_in_dir(2)
        x_coord = r_max + domain_min(1);
        for ii = 1:num_grains_in_dir(1)
            q(ctr,:) = [x_coord y_coord z_coord];
            x_coord = x_coord + d_max;
            ctr = ctr+1;
        end
        y_coord = y_coord + d_max;
    end
    z_coord = z_coord + d_max;
end

r = normrnd(r_mean,r_sigma,1,ctr-1);

end

