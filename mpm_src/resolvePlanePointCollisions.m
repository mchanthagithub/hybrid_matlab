function [ mpm_points ] = resolvePlanePointCollisions(mpm_points,static_planes)
tol = 1e-14;
[dummy,num_static_planes] = size(static_planes);
for pt_num = 1:mpm_points.num_points
    for plane_num = 1:num_static_planes
        dist = dot(static_planes(plane_num).n, mpm_points.q(pt_num,:) -static_planes(plane_num).q);
        if(dist > 0.0)
            continue;
        end
        
        % sliding
        if(static_planes(plane_num).type == 0)
            vnormal = dot(static_planes(plane_num).n,mpm_points.vel(pt_num,:));
            if(vnormal > 0.0)
                continue;
            end
            
            mpm_points.vel(pt_num,:) = mpm_points.vel(pt_num,:) - vnormal * static_planes(plane_num).n;
        % Sticking
        else
            mpm_points.vel(pt_num,:) = 0;
        end        
        
    end
    
end

end

