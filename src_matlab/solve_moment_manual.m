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
%   options: solver settings (cf. sdpsettings for YALMIP)
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
tic
NumMom = length(MomConst); NumLoc = length(LocConst); DimVar = size(obj, 2)-1;
constraints = [];
% check constant 1 variables
SetVars.var = sdpvar(size(obj, 1), 1); SetVars.supp = obj(:, 2:end);
turn = 'off';
for i = 1:length(SetVars.var)
    if isequal(obj(i, 2:end), zeros(1, DimVar))
%         sdisplay(SetVars.var(i))
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
% sdisplay(objective)
%
for i = 1:NumMom
    [SetVars, MomMat] = MomentMatrix(MomConst{i}.basis, MomConst{i}.ord, 'on', SetVars);
%     sdisplay(MomMat)
    constraints = [constraints, MomMat>=0];
    % update variable set
%     SetVars.var = [SetVars.var; vars.var];
%     SetVars.supp = [SetVars.supp; vars.supp];
    if isequal(turn, 'off')
        for k = 1:length(SetVars.var)
            if isequal(SetVars.supp(k, :), zeros(1, DimVar))
%                 sdisplay(SetVars.var(k))
                constraints = [constraints, SetVars.var(k)==1];
                turn = 'on';
                break
            end
        end
    end
%     for k = 1:length(vars.var)
%         for j = 1:length(SetVars.var)-length(vars.var)
%             if isequal(vars.supp(k,:), SetVars.supp(j,:))
%                 constraints = [constraints, vars.var(k)==SetVars.var(j)];
%                 break
%             end
%         end
%     end
end
%
for i = 1:NumLoc
    [SetVars, LocMat] = LocalizationMatrix(LocConst{i}.pol, LocConst{i}.basis, LocConst{i}.ord, 'on', SetVars);
%     sdisplay(LocMat)
    if isequal(LocConst{i}.typ, '>=')
        constraints = [constraints, LocMat>=0];
    elseif isequal(LocConst{i}.typ, '<=')
        constraints = [constraints, LocMat<=0];
    elseif isequal(LocConst{i}.typ, '==')
        constraints = [constraints, LocMat==0];
    end
    % update variable set
%     SetVars.var = [SetVars.var; vars.var];
%     SetVars.supp = [SetVars.supp; vars.supp];
    if isequal(turn, 'off')
        for k = 1:length(SetVars.var)
            if isequal(SetVars.supp(k, :), zeros(1, DimVar))
%                 sdisplay(SetVars.var(k))
                constraints = [constraints, SetVars.var(k)==1];
                turn = 'on';
                break
            end
        end
    end
%     for k = 1:length(vars.var)
%         for j = 1:length(SetVars.var)-length(vars.var)
%             if isequal(vars.supp(k,:), SetVars.supp(j,:))
%                 constraints = [constraints, vars.var(k)==SetVars.var(j)];
%                 break
%             end
%         end
%     end
end
% sdisplay(SetVars.var)
% getvariables(SetVars.var(1:3))
% constraints
time.model = toc;
%
fprintf('Problem type: %s\n', typ)
fprintf('%d moment matrices, %d localization matrices, %d variables, %d constraints\n', NumMom, NumLoc, length(SetVars.var), length(constraints))
fprintf('Solver started.............solving.............')
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
%
if isequal(stat, 'Successfully solved (MOSEK)')
    fprintf('The problem is successfully solved! Optimal value is: %f\n', OptVal)
else
    fprintf('The solver encountered some issues: %s\n', stat)
end
fprintf('Total running time is: %f (modeling time: %f, yalmip time: %f, solver time: %f)\n', time.model+time.yalm+time.solv, time.model, time.yalm, time.solv)
end