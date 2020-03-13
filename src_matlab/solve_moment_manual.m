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
%       options.duplicated: if 'on' set marginals to be equal inside
%           moments and localizations, outside otherwise
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
find_idx = find(ismember(obj(:, 2:end), zeros(1, DimVar), 'row') == 1, 1);
if ~isempty(find_idx)
    SetVars.var(find_idx) = 1;
    turn = 'on';
end
%
if max(sum(obj(:,2:end),2)) > 2
    ObjQuad = 'no';
    if isequal(typ, 'max')
        objective = - obj(:, 1)'*SetVars.var(1:size(obj, 1), 1);
    elseif isequal(typ, 'min')
        objective = obj(:, 1)'*SetVars.var(1:size(obj, 1), 1);
    end
else
    ObjQuad = 'yes';
    if isequal(options.duplicated, 'off')
        if isequal(typ, 'max')
            objective = - obj(:, 1)'*SetVars.var(1:size(obj, 1), 1);
        elseif isequal(typ, 'min')
            objective = obj(:, 1)'*SetVars.var(1:size(obj, 1), 1);
        end
    end
end
%
for i = 1:NumMom
    if i == 1
        fprintf('Building up moment matrices: %.2f%% (%d/%d)\n', i/NumMom*100, i, NumMom)
    else
        s1 = num2str(floor((i-1)/NumMom*100)); s2 = num2str(i-1); s3 = num2str(NumMom);
        fprintf([repmat('\b', 1, 9+length([s1,s2,s3])), '%.2f%% (%d/%d)\n'], i/NumMom*100, i, NumMom)
    end
    if isequal(options.duplicated, 'on')
        if i == 1
            [SetVars, MomMat] = MomentMatrix(MomConst{i}.basis, MomConst{i}.ord, options.duplicated, SetVars, ObjQuad);
            if isequal(typ, 'max')
                objective = - obj(:, 1)'*SetVars.var(1:size(obj, 1), 1);
            elseif isequal(typ, 'min')
                objective = obj(:, 1)'*SetVars.var(1:size(obj, 1), 1);
            end
        else
            [SetVars, MomMat] = MomentMatrix(MomConst{i}.basis, MomConst{i}.ord, options.duplicated, SetVars, 'no');
        end
%         sdisplay(MomMat)
        constraints = [constraints, MomMat>=0];
        if isequal(turn, 'off')
            find_idx = find(ismember(SetVars.supp, zeros(1, DimVar), 'row') == 1, 1);
            if ~isempty(find_idx)
                constraints = [constraints, SetVars.var(find_idx)==1];
                turn = 'on';
            end
        end
    elseif isequal(options.duplicated, 'off')
        [vars, MomMat] = MomentMatrix(MomConst{i}.basis, MomConst{i}.ord, options.duplicated, SetVars, 'no');
        constraints = [constraints, MomMat>=0];
        % update variable set
        SetVars.var = [SetVars.var; vars.var];
        SetVars.supp = [SetVars.supp; vars.supp];
        if isequal(turn, 'off')
            find_idx = find(ismember(SetVars.supp, zeros(1, DimVar), 'row') == 1, 1);
            if ~isempty(find_idx)
                constraints = [constraints, SetVars.var(find_idx)==1];
                turn = 'on';
            end
        end
        if isequal(ObjQuad, 'no') || (isequal(ObjQuad, 'yes') && length(vars.var) ~= length(SetVars.var)-length(vars.var))
            for k = 1:length(vars.var)
                find_idx = find(ismember(SetVars.supp(1:length(SetVars.var)-length(vars.var), :), vars.supp(k,:), 'row') == 1, 1);
                if ~isempty(find_idx)
                    constraints = [constraints, vars.var(k)==SetVars.var(find_idx)];
                    turn = 'on';
                end
            end
        elseif isequal(ObjQuad, 'yes') && i == 1 && length(vars.var) == length(SetVars.var)-length(vars.var)
            constraints = [constraints, vars.var==SetVars.var(1:length(vars.var))];
        end
    end
end
%
for i = 1:NumLoc
    if i == 1
        fprintf('Building up localization matrices: %.2f%% (%d/%d)\n', i/NumLoc*100, i, NumLoc)
    else
        s1 = num2str(floor((i-1)/NumLoc*100)); s2 = num2str(i-1); s3 = num2str(NumLoc);
        fprintf([repmat('\b', 1, 9+length([s1,s2,s3])), '%.2f%% (%d/%d)\n'], i/NumLoc*100, i, NumLoc)
    end
    if isequal(options.duplicated, 'on')
        [SetVars, LocMat] = LocalizationMatrix(LocConst{i}.pol, LocConst{i}.basis, LocConst{i}.ord, options.duplicated, SetVars);
%         sdisplay(LocMat)
        if isequal(LocConst{i}.typ, '>=')
            constraints = [constraints, LocMat>=0];
        elseif isequal(LocConst{i}.typ, '<=')
            constraints = [constraints, LocMat<=0];
        elseif isequal(LocConst{i}.typ, '==')
            constraints = [constraints, LocMat==0];
        end
        if isequal(turn, 'off')
            find_idx = find(ismember(SetVars.supp, zeros(1, DimVar), 'row') == 1, 1);
            if ~isempty(find_idx)
                constraints = [constraints, SetVars.var(find_idx)==1];
                turn = 'on';
            end
        end
    elseif isequal(options.duplicated, 'off')
        [vars, LocMat] = LocalizationMatrix(LocConst{i}.pol, LocConst{i}.basis, LocConst{i}.ord, options.duplicated, SetVars);
        if isequal(LocConst{i}.typ, '>=')
            constraints = [constraints, LocMat>=0];
        elseif isequal(LocConst{i}.typ, '<=')
            constraints = [constraints, LocMat<=0];
        elseif isequal(LocConst{i}.typ, '==')
            constraints = [constraints, LocMat==0];
        end
        % update variable set
        SetVars.var = [SetVars.var; vars.var];
        SetVars.supp = [SetVars.supp; vars.supp];
        if isequal(turn, 'off')
            find_idx = find(ismember(SetVars.supp, zeros(1, DimVar), 'row') == 1, 1);
            if ~isempty(find_idx)
                constraints = [constraints, SetVars.var(find_idx)==1];
                turn = 'on';
            end
        end
        if length(SetVars.var) ~= length(vars.var)
            for k = 1:length(vars.var)
                find_idx = find(ismember(SetVars.supp(1:length(SetVars.var)-length(vars.var), :), vars.supp(k,:), 'row') == 1, 1);
                if ~isempty(find_idx)
                    constraints = [constraints, vars.var(k)==SetVars.var(find_idx)];
                    turn = 'on';
                end
            end
        end
    end
end
% sdisplay(objective)
%
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
fprintf('Total running time is: %f seconds (modeling time: %f seconds, yalmip time: %f seconds, solver time: %f seconds)\n\n', time.model+time.yalm+time.solv, time.model, time.yalm, time.solv)
% double(SetVars.var)
end