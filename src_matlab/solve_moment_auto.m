function [OptVal, time, stat] = solve_moment_auto(typ, var, obj, MomConst, LocConst, options)
%% Solve general multi-order moment problem
%%
%
% INPUT:
%   typ: 'min' or 'max' (str)
%   var: basis of variables
%   obj: objective (sdpvar)
%   MomConst: moment sdp constraints (cell --> struct)
%       MomConst{i}.basis: monomial basis (sdpvar)
%       MomConst{i}.ord: order of moment matrix (double)
%   LocConst: localization sdp constraints (cell --> struct)
%       LocConst{i}.pol: constraints gi (sdpvar)
%       LocConst{i}.basis: monomial basis (sdpvar)
%       LocConst{i}.typ: '>=' or '<=' or '==' (str)
%       LocConst{i}.ord: order of localization matrix (double)
% 
% OUTPUT:
%   OptVal: optimaal value (double)
%   time: running time (struct)
%       time.solv: solving time
%       time.yalm: modeling time
%   stat: solution status
% 
%% Author: T. Chen
%%
obj = pol2supp(obj, var);
NumMom = length(MomConst); NumLoc = length(LocConst);
for i = 1:NumMom
    MomConst{i}.basis = basis2supp(MomConst{i}.basis, var);
end
for i = 1:NumLoc
    LocConst{i}.pol = pol2supp(LocConst{i}.pol, var);
    LocConst{i}.basis = basis2supp(LocConst{i}.basis, var);
end
[OptVal, time, stat] = solve_moment_manual(typ, obj, MomConst, LocConst, options);
end