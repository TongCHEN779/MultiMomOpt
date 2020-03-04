function [OptVal, time, stat] = solve_moment_manual(typ, obj, MomConst, LocConst, options)
%% Solve general multi-order moment problem
%%
%
% INPUT:
%   typ: 'min' or 'max' (str)
%   obj: [f_alpha | alpha], coefficients of objective w.r.t support 
%       (double)
%   MomConst: moment sdp constraints (cell --> struct)
%       MomConst{i}.basis: support of variables (double)
%       MomConst{i}.ord: order of moment matrix (double)
%   LocConst: localization sdp constraints (cell --> struct)
%       LocConst{i}.pol: [gi_alpha | alpha], coefficients of polynomial
%               w.r.t support (double)
%       LocConst{i}.basis: support of variables (double)
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
NumMom = length(MomConst); NumLoc = length(LocConst); DimVar = size(obj, 2)-1;
constraints = [];
% check constant 1 variables
SetVars.var = sdpvar(size(obj, 1), 1); SetVars.supp = obj(:, 2:end);
turn = 'off';
for i = 1:length(SetVars.var)
    if isequal(obj(i, 2:end), zeros(1, DimVar))
        SetVars.var(i) = 1;
        turn = 'on';
        break
    end
end
%
if isequal(typ, 'max')
    objective = - obj(:, 1)'*SetVars.var(1:size(obj, 1), 1);
elseif isequal(typ, 'min')
    objective = obj(:, 1)'*SetVars.var(1:size(obj, 1), 1);
end
%
for i = 1:NumMom
    [vars, MomMat] = moment(MomConst{i}.basis, MomConst{i}.ord);
    constraints = [constraints, MomMat>=0];
    % check common variables
    SetVars.var = [SetVars.var; vars.var];
    SetVars.supp = [SetVars.supp; vars.supp];
    if isequal(turn, 'off')
        for k = 1:length(vars.var)
            if isequal(vars.supp(k, :), zeros(1, DimVar))
                constraints = [constraints, vars.var(k)==1];
                turn = 'on';
                break
            end
        end
    end
    for k = 1:length(vars.var)
        for j = 1:length(SetVars.var)-length(vars.var)
            if isequal(vars.supp(k,:), SetVars.supp(j,:))
                constraints = [constraints, vars.var(k)==SetVars.var(j)];
                break
            end
        end
    end
end
%
for i = 1:NumLoc
    [vars, LocMat] = localization(LocConst{i}.pol, LocConst{i}.basis, LocConst{i}.ord);
%     constraints = [constraints, eval(sprintf('LocMat%s0', LocConst{i}.typ))];
    if isequal(LocConst{i}.typ, '>=')
        constraints = [constraints, LocMat>=0];
    elseif isequal(LocConst{i}.typ, '<=')
        constraints = [constraints, LocMat<=0];
    elseif isequal(LocConst{i}.typ, '==')
        constraints = [constraints, LocMat==0];
    end
    % check common variables
    SetVars.var = [SetVars.var; vars.var];
    SetVars.supp = [SetVars.supp; vars.supp];
    if isequal(turn, 'off')
        for k = 1:length(vars.var)
            if isequal(vars.supp(k, :), zeros(1, DimVar))
                constraints = [constraints, vars.var(k)==1];
                turn = 'on';
                break
            end
        end
    end
    for k = 1:length(vars.var)
        for j = 1:length(SetVars.var)-length(vars.var)
            if isequal(vars.supp(k,:), SetVars.supp(j,:))
                constraints = [constraints, vars.var(k)==SetVars.var(j)];
                break
            end
        end
    end
end
%
fprintf('Problem type: %s\n', typ)
fprintf('%d moment matrices, %d localization matrices, %d variables, %d constraints\n', NumMom, NumLoc, length(SetVars.var), length(constraints))
fprintf('Solver started...\n')
%
sol = optimize(constraints, objective, options);
fprintf('Solver finished.\n')
if isequal(typ, 'max')
    OptVal = value(-objective);
elseif isequal(typ, 'min')
    OptVal = value(objective);
end
time.solv = sol.solvertime;
time.yalm = sol.yalmiptime; 
stat = sol.info;
end