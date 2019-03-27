function [ contact_sphere_plane ] = findSpherePlaneCollisionsBruteForce( contact_sphere_plane,contact_sphere_plane_old,static_planes,q,r,v,del_t )
[n_dem,~] = size(q);
[n_static_planes,~] = size(static_planes);
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

end

