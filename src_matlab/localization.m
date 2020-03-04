function [vars, matrix] = localization(pol, basis, order)
%% Compute the generalized localization matix
%%
% 
% INPUT:
%   pol: [g_alpha | alpha], coefficients of polynomial w.r.t support 
%       (double)
%   basis: support of variables (double)
%   order: order of moment matrix (double)
%
% OUTPUT:
%   vars: variables needed in the localization matirx (struct)
%   matrix: localization matrix desired (sdpvar)
%
%% Author: T. Chen
%%
n = size(basis, 1);
B = get_basis(n,order); s = size(B,1); C = B*basis;
vars.var = []; vars.supp = sparse(0,0); matrix = 0;
de = size(pol, 1);
for i = 1:s
    for j = 1:s
        for k = 1:de
            new = 1;
            %
            for z = 1:length(vars.var)
                if isequal(vars.supp(z,:), pol(k,2:end)+C(i,:)+C(j,:))
                    E = zeros(s,s); E(i,j) = 1; E = sparse(E);
                    matrix = matrix + pol(k,1)*vars.var(z)*E;
                    new = 0;
                    break
                end
            end
            %
            if new == 1
                vars.var = [vars.var; sdpvar(1,1)];
                vars.supp = [vars.supp; sparse(pol(k,2:end)+C(i,:)+C(j,:))];
                E = zeros(s,s); E(i,j) = 1; E = sparse(E);
                matrix = matrix + pol(k,1)*vars.var(end)*E;
            end
        end
    end
end
end