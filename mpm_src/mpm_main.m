%% Code to 3D MPM
clc
clear
close all

%% Simulation parameters
t_final = 10.0;
del_t = 0.001;
L = 10;
grid_min = [0,0,0]; %x,y,z
grid_max = [L,L,L]; %x,y,z
grid_delta = L/20;
g = [0 0 -9.81];

%% Material properties
rho = 1000;
% Bulk modulus
K = 70e9;
% Poisson ratio
nu = 0.3;
% Shear modulus
G = 3*K*(1-2*nu)/(2*(1+nu));
points_per_dim = 2;

%% Make bodies
body(1).min = [5 5 5];
body(1).max = [7 7 7];

body(2).min = [2 5 5];
body(2).max = [4 7 7];

%% Generate Points
[dummy,num_bodies] = size(body);

num_points = 0;
point_delta = grid_delta / (points_per_dim+1);
for bdy_num = 1:num_bodies
    elem_min_idx = floor((body(bdy_num).min - grid_min)./grid_delta);
    elem_max_idx = ceil((body(bdy_num).max - grid_min)./grid_delta);
   
    num_points_prev = num_points;
    num_new_points = prod(elem_max_idx+1-elem_min_idx)*points_per_dim^3;
    num_points = num_points + num_new_points;
    
    q((num_points_prev+1):num_points,:) = zeros(num_new_points,3);
    
    coord_idx = num_points_prev+1;
    ctr = 0;
    for z_idx = elem_min_idx(3):elem_max_idx(3)
        for y_idx = elem_min_idx(2):elem_max_idx(2)
            for x_idx = elem_min_idx(1):elem_max_idx(1)
                base_coords = [x_idx-1 y_idx-1 z_idx-1]*grid_delta + grid_min;
                
                for z_inner_idx = 1:points_per_dim
                    for y_inner_idx = 1:points_per_dim
                        for x_inner_idx = 1:points_per_dim
                            q(coord_idx,:) = base_coords + point_delta*[x_inner_idx y_inner_idx z_inner_idx];
                            coord_idx = coord_idx+1;
                        end
                    end
                end
            end
        end
    end
end

scatter3(q(:,1),q(:,2),q(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
