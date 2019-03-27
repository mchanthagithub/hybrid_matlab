function [F,contact_sphere_plane] = calculateSpherePlaneForces(F,contact_sphere_plane,mass,k_n,k_t,gamma_n,gamma_t,mu)
    % Loop over sphere-plane contacts
    for kk = 1:contact_sphere_plane.num
        m_eff = mass(contact_sphere_plane.idx(kk,1));
        n = contact_sphere_plane.normal(kk,:);
        delta_s = contact_sphere_plane.delta_s(kk,:);
        pen_depth = contact_sphere_plane.pen_depth(kk);
        vel = contact_sphere_plane.v_rel(kk,:);
        vt = vel - dot(n,vel)*n;
        vn = vel - vt;
        
        normal_force = -k_n * pen_depth * n - ...
            0.5 * gamma_n * m_eff * vn;
        
        friction_force = -k_t * delta_s - 0.5 * m_eff * gamma_t * vt;
        
        % Project into friction cone
        mu_fn = mu * norm(normal_force);
        ft = norm(friction_force);
        if(0.5 * m_eff * gamma_t * norm(vt) > mu_fn)
           contact_sphere_plane.delta_s(kk,:) = 0;
           friction_force = -0.5 * m_eff * gamma_t * vt;
           friction_force = friction_force./norm(friction_force) * mu_fn;
        elseif(ft > mu_fn)
            a = k_t * k_t * norm(delta_s)^2;
            b = k_t * gamma_t * m_eff * dot(delta_s,vt);
            c = 0.25 * gamma_t * gamma_t * m_eff * m_eff * dot(vt,vt) - mu * mu * norm(normal_force)^2;
            dscr = b * b - 4.0 * a * c;
            dscr_sqrt = sqrt(dscr);
            root1 = secondRootOfQuadratic(a,b,c,dscr_sqrt);
            new_delta_s = root1 * delta_s;
            new_friction_force = -k_t * new_delta_s - 0.5 * m_eff * gamma_t * vt;
            friction_force = new_friction_force;
            contact_sphere_plane.delta_s(kk,:) = new_delta_s;
        end
        F(contact_sphere_plane.idx(kk),:) = F(contact_sphere_plane.idx(kk),:) + normal_force + friction_force;
    end

