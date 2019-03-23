%% 3D FEM that solves Poisson Equation (Assuming only Dirichlet and no neumann
% boundary conditions
% Assumes cube elements

clear
close all
clc

addpath('../mpm_src')

basis_functions = createBasisFunctions;

[connectivity, global_coordinates] = createElements(1,1,1,0.5,0.5,0.5)
