function [contact_sphere_sphere,contact_sphere_plane] = intializeContacts()
    contact_sphere_sphere.idx = zeros(0,0);
    contact_sphere_sphere.point = zeros(0,0);
    contact_sphere_sphere.normal = zeros(0,0);
    contact_sphere_sphere.pen_depth = zeros(0);
    contact_sphere_sphere.v_rel = zeros(0,0);
    contact_sphere_sphere.delta_s = zeros(0,0);
    contact_sphere_sphere.num = 0;
    
    contact_sphere_plane.idx = zeros(0,0);
    contact_sphere_plane.point = zeros(0,0);
    contact_sphere_plane.normal = zeros(0,0);
    contact_sphere_plane.pen_depth = zeros(0);
    contact_sphere_plane.v_rel = zeros(0,0);
    contact_sphere_plane.delta_s = zeros(0,0);
    contact_sphere_plane.num = 0;
end

