function [ basis_funcs ] = createBasisFunctions( )
% Only does linear basis functions for now

% Shape functions
basis_funcs.N{1} = @(xi) 0.125*(1-xi(1))*(1-xi(2))*(1-xi(3));
basis_funcs.N{2} = @(xi) 0.125*(1+xi(1))*(1-xi(2))*(1-xi(3));
basis_funcs.N{3} = @(xi) 0.125*(1+xi(1))*(1+xi(2))*(1-xi(3));
basis_funcs.N{4} = @(xi) 0.125*(1-xi(1))*(1+xi(2))*(1-xi(3));
basis_funcs.N{5} = @(xi) 0.125*(1-xi(1))*(1-xi(2))*(1+xi(3));
basis_funcs.N{6} = @(xi) 0.125*(1+xi(1))*(1-xi(2))*(1+xi(3));
basis_funcs.N{7} = @(xi) 0.125*(1+xi(1))*(1+xi(2))*(1+xi(3));
basis_funcs.N{8} = @(xi) 0.125*(1-xi(1))*(1+xi(2))*(1+xi(3));

% Gradients
% Bottom of element
basis_funcs.dN{1,1} = @(xi) -0.125*(1-xi(2))*(1-xi(3));
basis_funcs.dN{1,2} = @(xi) -0.125*(1-xi(1))*(1-xi(3));
basis_funcs.dN{1,3} = @(xi) -0.125*(1-xi(1))*(1-xi(2));

basis_funcs.dN{2,1} = @(xi) 0.125*(1-xi(2))*(1-xi(3));
basis_funcs.dN{2,2} = @(xi) -0.125*(1-xi(1))*(1-xi(3));
basis_funcs.dN{2,3} = @(xi) -0.125*(1-xi(1))*(1-xi(2));

basis_funcs.dN{3,1} = @(xi) 0.125*(1-xi(2))*(1-xi(3));
basis_funcs.dN{3,2} = @(xi) 0.125*(1-xi(1))*(1-xi(3));
basis_funcs.dN{3,3} = @(xi) -0.125*(1-xi(1))*(1-xi(2));

basis_funcs.dN{4,1} = @(xi) -0.125*(1-xi(2))*(1-xi(3));
basis_funcs.dN{4,2} = @(xi) 0.125*(1-xi(1))*(1-xi(3));
basis_funcs.dN{4,3} = @(xi) -0.125*(1-xi(1))*(1-xi(2));

% Top of element
basis_funcs.dN{5,1} = @(xi) -0.125*(1-xi(2))*(1-xi(3));
basis_funcs.dN{5,2} = @(xi) -0.125*(1-xi(1))*(1-xi(3));
basis_funcs.dN{5,3} = @(xi) 0.125*(1-xi(1))*(1-xi(2));

basis_funcs.dN{6,1} = @(xi) 0.125*(1-xi(2))*(1-xi(3));
basis_funcs.dN{6,2} = @(xi) -0.125*(1-xi(1))*(1-xi(3));
basis_funcs.dN{6,3} = @(xi) 0.125*(1-xi(1))*(1-xi(2));

basis_funcs.dN{7,1} = @(xi) 0.125*(1-xi(2))*(1-xi(3));
basis_funcs.dN{7,2} = @(xi) 0.125*(1-xi(1))*(1-xi(3));
basis_funcs.dN{7,3} = @(xi) 0.125*(1-xi(1))*(1-xi(2));

basis_funcs.dN{8,1} = @(xi) -0.125*(1-xi(2))*(1-xi(3));
basis_funcs.dN{8,2} = @(xi) 0.125*(1-xi(1))*(1-xi(3));
basis_funcs.dN{8,3} = @(xi) 0.125*(1-xi(1))*(1-xi(2));

end

