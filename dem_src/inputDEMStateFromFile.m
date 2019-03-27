function [q,r,v,material_parameters,contact_sphere_sphere,contact_sphere_plane,static_planes] = inputDEMStateFromFile(state_filename,contact_filename)
fprintf('Reading from %s and %s',state_filename,contact_filename);
state_fp = fopen(state_filename);
tline = fgetl(state_fp);
num_points = 0;
q_read_flag = 0;
v_read_flag = 0;
r_read_flag = 0;
material_param_read_flag = 0;
ctr = 0;
while ischar(tline)
%     disp(tline)
    
    if(q_read_flag == 1)
        ctr = ctr+1;
        tline = regexprep(tline,' +',' ');
        tline = strtrim(tline);
        tline = split(tline,' ');
        q(ctr,1) = str2num(tline{1});
        q(ctr,2) = str2num(tline{2});
        q(ctr,3) = str2num(tline{3});
        if(ctr == num_points)
            q_read_flag = 0;
            ctr = 0;
        end
    end
    
    if(v_read_flag == 1)
        ctr = ctr+1;
        tline = regexprep(tline,' +',' ');
        tline = strtrim(tline);
        tline = split(tline,' ');
        v(ctr,1) = str2num(tline{1});
        v(ctr,2) = str2num(tline{2});
        v(ctr,3) = str2num(tline{3});
        if(ctr == num_points)
            v_read_flag = 0;
            ctr = 0;
        end
    end
    
    if(r_read_flag == 1)
        ctr = ctr+1;
        tline = regexprep(tline,' +',' ');
        tline = strtrim(tline);
        tline = split(tline,' ');
        r(ctr,1) = str2num(tline{1});
        if(ctr == num_points)
            r_read_flag = 0;
            ctr = 0;
        end
    end
    
    if(material_param_read_flag == 1)
        ctr = ctr+1;
        tline = regexprep(tline,' +',' ');
        tline = strtrim(tline);
        tline = split(tline,' ');
        material_parameters(1) = str2num(tline{1});
        material_parameters(2) = str2num(tline{2});
        material_parameters(3) = str2num(tline{3});
        material_parameters(4) = str2num(tline{4});
        material_parameters(5) = str2num(tline{5});
        material_parameters(6) = str2num(tline{6});
        material_param_read_flag = 0;
    end
    
    if(contains(tline,'NumberOfPoints'))
        tline = split(tline,'"');
        num_points = str2num(tline{2});
        q = zeros(num_points,3);
        v = zeros(num_points,3);
        r = zeros(num_points,1);
    end
    
    if(contains(tline,'Position'))
        q_read_flag = 1;
    end
    
    if(contains(tline,'Velocity'))
        v_read_flag = 1;
    end
    
    if(contains(tline,'Radius'))
        r_read_flag = 1;
    end
    
    if(contains(tline,'Material_parameters'))
        material_param_read_flag = 1;
    end
    
    tline = fgetl(state_fp);
end
fclose(state_fp);


contact_fp = fopen(contact_filename);
tline = fgetl(contact_fp);
contact_sphere_sphere.num = 0;
contact_sphere_plane.num = 0;
idx_read_flag = 0;
point_read_flag = 0;
normal_read_flag = 0;
pen_depth_read_flag = 0;
v_rel_read_flag = 0;
delta_s_read_flag = 0;
sphere_sphere_flag = 0;
sphere_plane_flag = 0;
static_plane_flag = 0;
static_plane_q_flag = 0;
static_plane_n_flag = 0;
ctr = 0;
num = 0;
while ischar(tline)
%     disp(tline)
    
    if(idx_read_flag == 1)
        ctr = ctr+1;
        tline = split(tline,',');
        if(sphere_sphere_flag == 1)
           contact_sphere_sphere.idx(ctr,1) = str2num(tline{1});
           contact_sphere_sphere.idx(ctr,2) = str2num(tline{2});
        elseif(sphere_plane_flag == 1)
           contact_sphere_plane.idx(ctr,1) = str2num(tline{1});
           contact_sphere_plane.idx(ctr,2) = str2num(tline{2});
        end
        if(ctr == num)
            idx_read_flag = 0;
            ctr = 0;
        end
    end
    
    if(point_read_flag == 1)
        ctr = ctr+1;
        tline = split(tline,',');
        if(sphere_sphere_flag == 1)
           contact_sphere_sphere.point(ctr,1) = str2num(tline{1});
           contact_sphere_sphere.point(ctr,2) = str2num(tline{2});
           contact_sphere_sphere.point(ctr,3) = str2num(tline{3});
        elseif(sphere_plane_flag == 1)
           contact_sphere_plane.point(ctr,1) = str2num(tline{1});
           contact_sphere_plane.point(ctr,2) = str2num(tline{2});
           contact_sphere_plane.point(ctr,3) = str2num(tline{3});
        end
        if(ctr == num)
            point_read_flag = 0;
            ctr = 0;
        end
    end
    
    if(normal_read_flag == 1)
        ctr = ctr+1;
        tline = split(tline,',');
        if(sphere_sphere_flag == 1)
           contact_sphere_sphere.normal(ctr,1) = str2num(tline{1});
           contact_sphere_sphere.normal(ctr,2) = str2num(tline{2});
           contact_sphere_sphere.normal(ctr,3) = str2num(tline{3});
        elseif(sphere_plane_flag == 1)
           contact_sphere_plane.normal(ctr,1) = str2num(tline{1});
           contact_sphere_plane.normal(ctr,2) = str2num(tline{2});
           contact_sphere_plane.normal(ctr,3) = str2num(tline{3});
        end
        if(ctr == num)
            normal_read_flag = 0;
            ctr = 0;
        end
    end
    
    if(pen_depth_read_flag == 1)
        ctr = ctr+1;
        if(sphere_sphere_flag == 1)
           contact_sphere_sphere.pen_depth(ctr,1) = str2num(tline);
        elseif(sphere_plane_flag == 1)
           contact_sphere_plane.pen_depth(ctr,1) = str2num(tline);
        end
        if(ctr == num)
            pen_depth_read_flag = 0;
            ctr = 0;
        end
    end
    
    if(v_rel_read_flag == 1)
        ctr = ctr+1;
        tline = split(tline,',');
        if(sphere_sphere_flag == 1)
           contact_sphere_sphere.v_rel(ctr,1) = str2num(tline{1});
           contact_sphere_sphere.v_rel(ctr,2) = str2num(tline{2});
           contact_sphere_sphere.v_rel(ctr,3) = str2num(tline{3});
        elseif(sphere_plane_flag == 1)
           contact_sphere_plane.v_rel(ctr,1) = str2num(tline{1});
           contact_sphere_plane.v_rel(ctr,2) = str2num(tline{2});
           contact_sphere_plane.v_rel(ctr,3) = str2num(tline{3});
        end
        if(ctr == num)
            v_rel_read_flag = 0;
            ctr = 0;
        end
    end
    
    if(delta_s_read_flag == 1)
        ctr = ctr+1;
        tline = split(tline,',');
        if(sphere_sphere_flag == 1)
           contact_sphere_sphere.delta_s(ctr,1) = str2num(tline{1});
           contact_sphere_sphere.delta_s(ctr,2) = str2num(tline{2});
           contact_sphere_sphere.delta_s(ctr,3) = str2num(tline{3});
        elseif(sphere_plane_flag == 1)
           contact_sphere_plane.delta_s(ctr,1) = str2num(tline{1});
           contact_sphere_plane.delta_s(ctr,2) = str2num(tline{2});
           contact_sphere_plane.delta_s(ctr,3) = str2num(tline{3});
        end
        if(ctr == num)
            delta_s_read_flag = 0;
            ctr = 0;
        end
    end
    
    if(static_plane_q_flag == 1)
        ctr = ctr+1;
        tline = split(tline,',');
        if(static_plane_flag == 1)
           temp_q = [str2num(tline{1}) str2num(tline{2}) str2num(tline{3})];
           static_planes(ctr).q = temp_q;
        end
        if(ctr == num)
            static_plane_q_flag = 0;
            ctr = 0;
        end
    end
    
    
    if(static_plane_n_flag == 1)
        ctr = ctr+1;
        tline = split(tline,',');
        if(static_plane_flag == 1)
           temp_n = [str2num(tline{1}) str2num(tline{2}) str2num(tline{3})];
           static_planes(ctr).n = temp_n;
        end
        if(ctr == num)
            static_plane_n_flag = 0;
            ctr = 0;
        end
    end
    
    if(contains(tline,'num'))
        tline_old = tline;
        tline_old = split(tline_old,' ,');
        tline = fgetl(contact_fp);
        num = str2num(tline);
        if(contains(tline_old{1},'sphere_sphere'))
            sphere_sphere_flag = 1;
            contact_sphere_sphere.num = num;
            contact_sphere_sphere.idx = zeros(num,2);
            contact_sphere_sphere.point = zeros(num,3);
            contact_sphere_sphere.normal = zeros(num,3);
            contact_sphere_sphere.pen_depth = zeros(num,1);
            contact_sphere_sphere.v_rel = zeros(num,3);
            contact_sphere_sphere.delta_s = zeros(num,3);
            if(num == 0)
                sphere_sphere_flag = 0;
            end
        elseif(contains(tline_old{1},'sphere_plane'))
            sphere_plane_flag = 1;
            sphere_sphere_flag = 0;
            contact_sphere_plane.num = num;
            contact_sphere_plane.idx = zeros(num,2);
            contact_sphere_plane.point = zeros(num,3);
            contact_sphere_plane.normal = zeros(num,3);
            contact_sphere_plane.pen_depth = zeros(num,1);
            contact_sphere_plane.v_rel = zeros(num,3);
            contact_sphere_plane.delta_s = zeros(num,3);
            if(num == 0)
                sphere_plane_flag = 0;
            end
        elseif(contains(tline_old{1},'static_planes'))
            sphere_plane_flag = 0;
            sphere_sphere_flag = 0;
            static_plane_flag = 1;
            static_planes.num = num;
            static_planes.q = zeros(num,3);
            static_planes.n = zeros(num,3);
            if(num == 0)
                static_plane_flag = 0;
            end
        end
    end
    

    if((sphere_sphere_flag || sphere_plane_flag) && sum(contains(tline,'idx')))
        idx_read_flag = 1;
    end
    
    if((sphere_sphere_flag || sphere_plane_flag) && sum(contains(tline,'point')))
        point_read_flag = 1;
    end
    
    if((sphere_sphere_flag || sphere_plane_flag) && sum(contains(tline,'normal')))
        normal_read_flag = 1;
    end
    
    if((sphere_sphere_flag || sphere_plane_flag) && sum(contains(tline,'pen_depth')))
        pen_depth_read_flag = 1;
    end
    
    if((sphere_sphere_flag || sphere_plane_flag) && sum(contains(tline,'v_rel')))
        v_rel_read_flag = 1;
    end
    
    if((sphere_sphere_flag || sphere_plane_flag) && sum(contains(tline,'delta_s')))
        delta_s_read_flag = 1;
    end
    
    if(static_plane_flag && sum(contains(tline,'q')))
        static_plane_q_flag = 1;
    end
    
    if(static_plane_flag && sum(contains(tline,'normal')))
        static_plane_n_flag = 1;
    end
    
    tline = fgetl(contact_fp);
end
fclose(contact_fp);
