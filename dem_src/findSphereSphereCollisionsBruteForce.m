function [ contact_sphere_sphere ] = findSphereSphereCollisionsBruteForce( contact_sphere_sphere,contact_sphere_sphere_old,q,r,v,del_t )
[n_dem, dummy] = size(q);
for kk = 1:n_dem
        for jj = kk+1:n_dem
            diff = q(kk,:) - q(jj,:);
            dist = norm(diff);
            if(dist < r(jj) + r(kk))
                contact_sphere_sphere.num = contact_sphere_sphere.num + 1;
                idx = [kk jj];
                normal = diff./dist;
                pen_depth = dist - (r(jj) + r(kk));
                point = q(kk,:) + normal*pen_depth/2.0;
                v_rel = v(kk,:) - v(jj,:);
                contact_sphere_sphere.idx(contact_sphere_sphere.num,:) = idx;
                contact_sphere_sphere.normal(contact_sphere_sphere.num,:) = normal;
                contact_sphere_sphere.pen_depth(contact_sphere_sphere.num) = pen_depth;
                contact_sphere_sphere.point(contact_sphere_sphere.num,:) = point;
                contact_sphere_sphere.v_rel(contact_sphere_sphere.num,:) = v_rel;
                contact_exists = 0;
                if(contact_sphere_sphere_old.num > 0)
                    [contact_exists, index] = ismember([kk jj],contact_sphere_sphere_old.idx,'rows');
                end
                if(contact_exists)
                    delta_s = contact_sphere_sphere_old.delta_s(index,:);
                    delta_s = delta_s + del_t*v_rel;
                    delta_s = delta_s - dot(normal,delta_s)*normal;
                    contact_sphere_sphere.delta_s(contact_sphere_sphere.num,:) = delta_s;
                else
                    delta_s = del_t * v_rel;
                    delta_s = delta_s - dot(normal,delta_s)*normal;
                    contact_sphere_sphere.delta_s(contact_sphere_sphere.num,:) = delta_s;
                end
            end
        end
    end


end

