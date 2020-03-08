function [vars, matrix] = MomentMatrix(basis, order, duplicated, SetVars, ObjQuad)
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
if isequal(duplicated, 'off')
    vars.var = []; vars.supp = [];
    if order == 1
        n = size(basis, 1);
        E = sparse([zeros(1+n, n-1), [2, zeros(1,n); ones(n,1), eye(n)]]);
        matrix = sdpvar(n+1);
        for i = 1:n+1
            vars.var = [vars.var; matrix(i:end, i)];
            vars.supp = [vars.supp; E(1:end-i+1, end-n-i+2:end-i+1)];
        end
    elseif order >= 2
        pol = [1, zeros(1, size(basis,2))];
        [vars, matrix] = LocalizationMatrix(pol, basis, order, duplicated, SetVars);
    end
elseif isequal(duplicated, 'on')
    if order == 1
        vars = SetVars;
        n = size(basis, 1);
        E = sparse([zeros(1+n, n-1), [2, zeros(1,n); ones(n,1), eye(n)]]);
        matrix = sdpvar(n+1);
        for i = 1:n+1
            % check common variables
            if isequal(ObjQuad, 'no')
                for j = i:n+1
                    new = 1;
                    find_idx = find(ismember(vars.supp, E(j-i+1, end-n-i+2:end-i+1), 'row') == 1, 1);
                    if ~isempty(find_idx)
                        matrix(j,i) = vars.var(find_idx); matrix(i,j) = vars.var(find_idx);
                        new = 0;
                    end
%                     for k = 1:length(vars.var)
%                         if isequal(E(j-i+1, end-n-i+2:end-i+1), vars.supp(k,:))
%                             matrix(j,i) = vars.var(k); matrix(i,j) = vars.var(k);
%                             new = 0;
%                             break
%                         end
%                     end
                    %
                    if new == 1
                        vars.var = [vars.var; matrix(j,i)];
                        vars.supp = [vars.supp; E(j-i+1, end-n-i+2:end-i+1)];
                    end
                end
            elseif isequal(ObjQuad, 'yes')
%                 i
                matrix(i:end, i) = vars.var((2*n+4-i)*(i-1)/2 + 1:(2*n+3-i)*(i)/2); 
                matrix(i, i:end) = vars.var((2*n+4-i)*(i-1)/2 + 1:(2*n+3-i)*(i)/2);
            end
        end
    elseif order >= 2
        pol = [1, zeros(1, size(basis,2))];
        [vars, matrix] = LocalizationMatrix(pol, basis, order, duplicated, SetVars);
    end
end
% n = length(basis.var);
% B = get_basis(n,2*order); s = size(B,1);
% newvars.var = sdpvar(s-n-1,1);
% newvars.supp = B(n+2:end,:)*basis.supp;
% variable = sprintf("moment_%d_%d", n, order);
% load('Moment_and_Localization_Matrices.mat',variable);
% Matrix_temp = eval(variable); Psd = [1; basis.var; newvars.var];
% Matrix = 0;
% for t = 1:s
%     Matrix = Matrix + Matrix_temp{t}*Psd(t);
% end
end