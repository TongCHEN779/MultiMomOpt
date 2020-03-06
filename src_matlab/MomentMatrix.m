function [vars, matrix] = MomentMatrix(basis, order, duplicated, SetVars)
%% Compute the moment matrix
%%
%
% INPUT:
%   basis: support of variables (double)
%   order: order of moment matrix (double)
%   duplicated: 'off' if duplicated variables are allowed, 'on' otherwise
%   SetVars: variable set of the outer environment (only makes sense when 
%       duplicated = 'on')
%
% OUTPUT:
%   vars: variables needed in the moment matirx (struct)
%   matrix: moment matrix desired (sdpvar)
%
%% Author: T. Chen
%%
if (nargin == 2) || (nargin == 4 && isequal(duplicated, 'off'))
    if order == 1
        n = size(basis, 1);
        E = sparse([zeros(1+n, n-1), [2, zeros(1,n); ones(n,1), eye(n)]]);
        matrix = sdpvar(n+1);
        vars.var = []; vars.supp = [];
        for i = 1:n+1
            vars.var = [vars.var; matrix(i:end, i)];
            vars.supp = [vars.supp; E(1:end-i+1, end-n-i+2:end-i+1)];
        end
    elseif order >= 2
        pol = [1, zeros(1, size(basis,2))];
        [vars, matrix] = localization(pol, basis, order);
    end
elseif nargin == 4 && isequal(duplicated, 'on')
    vars = SetVars;
    if order == 1
        n = size(basis, 1);
        E = sparse([zeros(1+n, n-1), [2, zeros(1,n); ones(n,1), eye(n)]]);
        matrix = sdpvar(n+1);
        for i = 1:n+1
            % check common variables
            for j = i:n+1
                new = 1;
                for k = 1:length(vars.var)
                    if isequal(E(j-i+1, end-n-i+2:end-i+1), vars.supp(k,:))
                        matrix(j,i) = vars.var(k); matrix(i,j) = vars.var(k);
                        new = 0;
                        break
                    end
                end
                %
                if new == 1
                    vars.var = [vars.var; matrix(j,i)];
                    vars.supp = [vars.supp; E(j-i+1, end-n-i+2:end-i+1)];
                end
            end
        end
    elseif order >= 2
        pol = [1, zeros(1, size(basis,2))];
        [vars, matrix] = LocalizationMatrix(pol, basis, order);
    end
end
end