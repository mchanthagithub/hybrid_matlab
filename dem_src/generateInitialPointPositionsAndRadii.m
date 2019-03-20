function [ q,r ] = generateInitialPointPositionsAndRadii( n_dem,r_mean,r_sigma,L,domain_min,domain_max,rngseed )
rng(rngseed)
q = zeros(0,0);
for kk = 1:n_dem
   overlap = true;
   [num_valid_grains,dim] = size(q);
   while(overlap)
       q_temp = L*rand(1,3);
       r_temp = normrnd(r_mean,r_sigma);
       
       
       
       overlap_temp = false;
       for ii = 1:num_valid_grains
          dist = sqrt( sum((q_temp-q(ii,:)).^2) );
          if(dist < r_temp + r(ii))
              overlap_temp = true;
              break;
          end
       end
       
       if(q_temp(1) - r_temp < domain_min(1) || q_temp(2) - r_temp < domain_min(2) || q_temp(3) - r_temp < domain_min(3))
           overlap_temp = true;
       end
       
        if(q_temp(1) + r_temp > domain_max(1) || q_temp(2) + r_temp > domain_max(2) || q_temp(3) + r_temp > domain_max(3))
           overlap_temp = true;
       end
       
       overlap = overlap_temp;
   end
   q(kk,:) = q_temp;
   r(kk) = r_temp;

end

end

