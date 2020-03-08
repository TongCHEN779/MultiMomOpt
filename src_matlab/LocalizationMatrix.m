function [vars, matrix] = LocalizationMatrix(pol, basis, order, duplicated, SetVars)
%% Compute the generalized localization matix
%%
%
% INPUT:
%   pol: [g_alpha | alpha], coefficients of polynomial w.r.t support
%       (double)
%   basis: support of variables (double)
%   order: order of moment matrix (double)
%   duplicated: 'off' if duplicated variables are allowed, 'on' otherwise
%   SetVars: variable set of the outer environment (only makes sense when 
%       duplicated = 'on')
%
% OUTPUT:
%   vars: variables needed in the localization matirx (struct)
%   matrix: localization matrix desired (sdpvar)
%
%% Author: T. Chen
%%
n = size(basis, 1);
B = get_basis(n,order); s = size(B,1); C = B*basis;
matrix = 0;
de = size(pol, 1);
if isequal(duplicated, 'off')
    vars.var = []; vars.supp = [];
    for i = 1:s
        for j = i:s
            for k = 1:de
                new = 1;
                %
                if ~isempty(vars.supp)
                    find_idx = find(ismember(vars.supp, pol(k,2:end)+C(i,:)+C(j,:), 'row') == 1, 1);
                    if ~isempty(find_idx)
                        E = zeros(s,s); E(i,j) = 1; E(j,i) = 1; E = sparse(E);
                        matrix = matrix + pol(k,1)*vars.var(find_idx)*E;
                        new = 0;
                    end
                end
                %
                if new == 1
                    vars.var = [vars.var; sdpvar(1,1)];
                    vars.supp = [vars.supp; sparse(pol(k,2:end)+C(i,:)+C(j,:))];
                    E = zeros(s,s); E(i,j) = 1; E(j,i) = 1; E = sparse(E);
                    matrix = matrix + pol(k,1)*vars.var(end)*E;
                end
            end
        end
    end
elseif isequal(duplicated, 'on')
    vars = SetVars;
    for i = 1:s
        for j = i:s
            for k = 1:de
                new = 1;
                % check common variables
                find_idx = find(ismember(vars.supp, pol(k,2:end)+C(i,:)+C(j,:), 'row') == 1, 1);
                if ~isempty(find_idx)
                    E = zeros(s,s); E(i,j) = 1; E(j,i) = 1; E = sparse(E);
                    matrix = matrix + pol(k,1)*vars.var(find_idx)*E;
                    new = 0;
                end
                %
                if new == 1
                    vars.var = [vars.var; sdpvar(1,1)];
                    vars.supp = [vars.supp; sparse(pol(k,2:end)+C(i,:)+C(j,:))];
                    E = zeros(s,s); E(i,j) = 1; E(j,i) = 1; E = sparse(E);
                    matrix = matrix + pol(k,1)*vars.var(end)*E;
                end
            end
        end
    end
end
end