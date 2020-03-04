function supp = pol2supp(pol, x)
%% Convert sdpvar polynomial object to double support object
%%
%
% INPUT:
%   pol: polynomial (sdpvar)
%   x: n-dimensional variable
% 
% OUTPUT:
%   supp: [f_alpha | alpha] coefficients and support of f (double)
%
%% Author: T. Chen
%%
[c, v] = coefficients(pol, x);
NumTerms = length(v); supp = sparse(0,0);
for i = 1:NumTerms
    supp = [supp; [c(i), degree(v(i), x)]];
end
end