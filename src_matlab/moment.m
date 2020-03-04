function [vars, matrix] = moment(basis, order)
%% Compute the moment matrix
%%
% 
% INPUT:
%   basis: support of variables (double)
%   order: order of moment matrix (double)
%
% OUTPUT:
%   vars: variables needed in the moment matirx (struct)
%   matrix: moment matrix desired (sdpvar)
%
%% Author: T. Chen
%%
pol = [1, zeros(1, size(basis,2))];
[vars, matrix] = localization(pol, basis, order);
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