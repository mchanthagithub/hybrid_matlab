function [ contact_sphere_sphere ] = findSphereSphereCollisionsUsingGrid( contact_sphere_sphere,contact_sphere_sphere_old, grid,q,r,v,del_t )
[num_grid_elements,dummy] = size(grid.collisions);
for elem_num = 1:num_grid_elements
    num_grains_in_cell = length(grid.collisions{elem_num});
    for kk = 1:num_grains_in_cell
        for jj = kk+1:num_grains_in_cell
            idx1 = grid.collisions{elem_num}(kk);
            idx2 = grid.collisions{elem_num}(jj);
%             contact_exists = 0;
%             if(contact_sphere_sphere.num > 0)
%                 [contact_exists, index] = ismember([idx1 idx2],contact_sphere_sphere.idx,'rows');
%             end
%             if(contact_exists)
%                 continue
%             end
            
            diff = q(idx1,:) - q(idx2,:);
            dist = norm(diff);
            if(dist < r(idx2) + r(idx1))
                contact_sphere_sphere.num = contact_sphere_sphere.num + 1;
                idx = [idx1 idx2];
                normal = diff./dist;
                pen_depth = dist - (r(idx2) + r(idx1));
                point = q(idx1,:) + normal*pen_depth/2.0;
                v_rel = v(idx1,:) - v(idx2,:);
                contact_sphere_sphere.idx(contact_sphere_sphere.num,:) = idx;
                contact_sphere_sphere.normal(contact_sphere_sphere.num,:) = normal;
                contact_sphere_sphere.pen_depth(contact_sphere_sphere.num) = pen_depth;
                contact_sphere_sphere.point(contact_sphere_sphere.num,:) = point;
                contact_sphere_sphere.v_rel(contact_sphere_sphere.num,:) = v_rel;
                contact_exists = 0;
                if(contact_sphere_sphere_old.num > 0)
                    [contact_exists, index] = ismember([idx1 idx2],contact_sphere_sphere_old.idx,'rows');
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

% Get rid of duplicate entries
[idx_new,map,~] = unique(contact_sphere_sphere.idx,'rows');
[new_contact_num, ~] = size(idx_new);
if(new_contact_num ~= contact_sphere_sphere.num)
    contact_sphere_sphere.num = new_contact_num;
    contact_sphere_sphere.idx = idx_new;
    contact_sphere_sphere.normal = contact_sphere_sphere.normal(map,:);
    contact_sphere_sphere.pen_depth = contact_sphere_sphere.pen_depth(map);
    contact_sphere_sphere.point = contact_sphere_sphere.point(map,:);
    contact_sphere_sphere.v_rel = contact_sphere_sphere.v_rel(map,:);
    contact_sphere_sphere.delta_s = contact_sphere_sphere.delta_s(map,:);
end



end

