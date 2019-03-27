function [] = outputDEMContacts(filename,contact_sphere_sphere,contact_sphere_plane,static_planes)

fp = fopen(filename,'w');

fprintf(fp,'contact_sphere_sphere, num\n');
fprintf(fp,'%d\n',contact_sphere_sphere.num);
fprintf(fp,'contact_sphere_sphere, idx\n');
for ii = 1:contact_sphere_sphere.num
   fprintf(fp,'%d,%d\n',contact_sphere_sphere.idx(ii,1),contact_sphere_sphere.idx(ii,2)); 
end

fprintf(fp,'contact_sphere_sphere, point\n');
for ii = 1:contact_sphere_sphere.num
   fprintf(fp,'%.16g,%.16g,%.16g\n',contact_sphere_sphere.point(ii,1),contact_sphere_sphere.point(ii,2),contact_sphere_sphere.point(ii,3)); 
end

fprintf(fp,'contact_sphere_sphere, normal\n');
for ii = 1:contact_sphere_sphere.num
   fprintf(fp,'%.16g,%.16g,%.16g\n',contact_sphere_sphere.normal(ii,1),contact_sphere_sphere.normal(ii,2),contact_sphere_sphere.normal(ii,3)); 
end

fprintf(fp,'contact_sphere_sphere, pen_depth\n');
for ii = 1:contact_sphere_sphere.num
   fprintf(fp,'%.16g\n',contact_sphere_sphere.pen_depth(ii)); 
end

fprintf(fp,'contact_sphere_sphere, v_rel\n');
for ii = 1:contact_sphere_sphere.num
   fprintf(fp,'%.16g,%.16g,%.16g\n',contact_sphere_sphere.v_rel(ii,1),contact_sphere_sphere.v_rel(ii,2),contact_sphere_sphere.v_rel(ii,3)); 
end

fprintf(fp,'contact_sphere_sphere, delta_s\n');
for ii = 1:contact_sphere_sphere.num
   fprintf(fp,'%.16g,%.16g,%.16g\n',contact_sphere_sphere.delta_s(ii,1),contact_sphere_sphere.delta_s(ii,2),contact_sphere_sphere.delta_s(ii,3)); 
end

%%

fprintf(fp,'contact_sphere_plane, num\n');
fprintf(fp,'%d\n',contact_sphere_plane.num);
fprintf(fp,'contact_sphere_plane, idx\n');
for ii = 1:contact_sphere_plane.num    
   fprintf(fp,'%d,%d\n',contact_sphere_plane.idx(ii,1),contact_sphere_plane.idx(ii,2)); 
end

fprintf(fp,'contact_sphere_plane, point\n');
for ii = 1:contact_sphere_plane.num
   fprintf(fp,'%.16g,%.16g,%.16g\n',contact_sphere_plane.point(ii,1),contact_sphere_plane.point(ii,2),contact_sphere_plane.point(ii,3)); 
end

fprintf(fp,'contact_sphere_plane, normal\n');
for ii = 1:contact_sphere_plane.num
   fprintf(fp,'%.16g,%.16g,%.16g\n',contact_sphere_plane.normal(ii,1),contact_sphere_plane.normal(ii,2),contact_sphere_plane.normal(ii,3)); 
end

fprintf(fp,'contact_sphere_plane, pen_depth\n');
for ii = 1:contact_sphere_plane.num
   fprintf(fp,'%.16g\n',contact_sphere_plane.pen_depth(ii)); 
end

fprintf(fp,'contact_sphere_plane, v_rel\n');
for ii = 1:contact_sphere_plane.num
   fprintf(fp,'%.16g,%.16g,%.16g\n',contact_sphere_plane.v_rel(ii,1),contact_sphere_plane.v_rel(ii,2),contact_sphere_plane.v_rel(ii,3)); 
end

fprintf(fp,'contact_sphere_plane, delta_s\n');
for ii = 1:contact_sphere_plane.num
   fprintf(fp,'%.16g,%.16g,%.16g\n',contact_sphere_plane.delta_s(ii,1),contact_sphere_plane.delta_s(ii,2),contact_sphere_plane.delta_s(ii,3)); 
end

[~,num_planes] = size(static_planes);
fprintf(fp,'static_planes, num\n');
fprintf(fp,'%d\n',num_planes);

fprintf(fp,'static_planes, q\n');
for ii = 1:num_planes
   fprintf(fp,'%.16g,%.16g,%.16g\n',static_planes(ii).q); 
end

fprintf(fp,'static_planes, normal\n');
for ii = 1:num_planes
   fprintf(fp,'%.16g,%.16g,%.16g\n',static_planes(ii).n); 
end


fclose(fp);


end

