function supp = basis2supp(basis, x)
%% Convert sdpvar monomial basis object to double support object
%%
%
% INPUT:
%   basis: monomial basis (sdpvar)
%   x: n-dimensional variable
% 
% OUTPUT:
%   supp: [f_alpha | alpha] coefficients and support of f (double)
%
%% Author: T. Chen
%%
NumBasis = length(basis); supp = sparse(0, 0);
for i = 1:NumBasis
    supp = [supp; degree(basis(i), x)];
end
end