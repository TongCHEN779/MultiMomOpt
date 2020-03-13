function [OptVal, time, stat, obj, MomConst, LocConst] = solve_moment_mip(A, options)
%% Solve MAX-CUT problem using rounding technique
%%
%
% INPUT:
%   A: weight matrix (double)
%   options: solver settings (cf. sdpsettings for YALMIP)
%       options.duplicated: if 'on' set marginals to be equal inside
%           moments and localizations, outside otherwise (str)
%       options.clique: if 'on' use subhierarchy w.r.t. cliques in the
%           chordal graph, use canonical order if 'off' (str)
%       options.order: relaxation order (double)
%       options.level: level of subhierarchy (double)
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
if nargin == 2
end
n = size(A, 1); typ = 'min';
%
obj = []; E0 = sparse([zeros(1+n, n-1), [2, zeros(1,n); ones(n,1), eye(n)]]); E = sparse([eye(n,n); eye(n,n)]);
A = 2*A;
for i = 1:n
    A(i,i) = A(i,i)/2;
end
obj = [obj, [zeros(n+1,1), E0(1:end, end-n+1:end)]];
for i = 2:n+1
    obj = [obj; [A(i-1:end,i-1), E0(1:end-i+1, end-n-i+2:end-i+1)]];
end
%
if options.level == 0 % Shor's relaxation
    MomConst{1}.basis = sparse(eye(n,n)); MomConst{1}.ord = 1;
    %
    for i = 1:n
        LocConst{i}.pol = sparse([1, E(i,:); -1, 2*E(i,:)]);
        LocConst{i}.basis = sparse(eye(n,n));
        LocConst{i}.typ = '==';
        LocConst{i}.ord = 0;
    end
elseif options.level > 0 && isequal(options.clique, 'on')
    [cliques, sizes] = chordal_extension(A, 'off');
    %
    MomConst{1}.basis = sparse(eye(n,n)); MomConst{1}.ord = 1;
    NumMom = 1;
    for i = 1:length(cliques)
        clique = [cliques{i}, cliques{i}];
        if options.level < sizes(i)
            for j = 1:sizes(i)
                NumMom = NumMom + 1;
                MomConst{NumMom}.basis = sparse(E(clique(j:j+options.level-1), :));
                MomConst{NumMom}.ord = options.order;
            end
        else
            NumMom = NumMom + 1;
            MomConst{NumMom}.basis = sparse(E(cliques{i}, :));
            MomConst{NumMom}.ord = options.order;
        end
    end
    %
    idx = [];
    NumLoc = 0;
    for i = 1:n
        for j = 1:length(cliques)
            if ismember(i, cliques{j})
                clique = [cliques{j}, cliques{j}]; idx = j;
                break
            end
        end
        if options.level < sizes(idx)
            for j = 1:sizes(idx)
                NumLoc = NumLoc + 1;
                LocConst{NumLoc}.pol = sparse([1, E(i,:); -1, 2*E(i,:)]);
                LocConst{NumLoc}.basis = sparse(E(clique(j:j+options.level-1), :));
                LocConst{NumLoc}.typ = '==';
                LocConst{NumLoc}.ord = options.order - 1;
            end
        else
            NumLoc = NumLoc + 1;
            LocConst{NumLoc}.pol = sparse([1, E(i,:); -1, 2*E(i,:)]);
            LocConst{NumLoc}.basis = sparse(E(cliques{idx}, :));
            LocConst{NumLoc}.typ = '==';
            LocConst{NumLoc}.ord = options.order - 1;
        end
    end
elseif options.level > 0 && isequal(options.clique, 'off')
    MomConst{1}.basis = sparse(eye(n,n));
    MomConst{1}.ord = 1;
    NumMom = 1; NumLoc = 0;
    clique = [1:n, 1:n];
    if options.level < n
        for i = 1:n
            NumMom = NumMom + 1;
            MomConst{NumMom}.basis = sparse(E(clique(i:i+options.level-1), :));
            MomConst{NumMom}.ord = options.order;
            %
            NumLoc = NumLoc + 1;
            LocConst{NumLoc}.pol = sparse([1, E(i,:); -1, 2*E(i,:)]);
            LocConst{NumLoc}.basis = sparse(E(clique(i:i+options.level-1), :));
            LocConst{NumLoc}.typ = '==';
            LocConst{NumLoc}.ord = options.order - 1;
        end
    else
        NumMom = NumMom + 1;
        MomConst{NumMom}.basis = sparse(E(1:n, :));
        MomConst{NumMom}.ord = options.order;
        %
        for i = 1:n
            NumLoc = NumLoc + 1;
            LocConst{NumLoc}.pol = sparse([1, E(i,:); -1, 2*E(i,:)]);
            LocConst{NumLoc}.basis = sparse(E(1:n, :));
            LocConst{NumLoc}.typ = '==';
            LocConst{NumLoc}.ord = options.order - 1;
        end
    end
end
%
fprintf('\nMixed Integer Programming (MIP): %d variables, %d entries, duplicated %s, clique %s, order %d, level %d\n', n, (sum(sum(full(1*(A~=0))))-sum(1*diag(A~=0)))/2+sum(1*diag(A~=0)), options.duplicated, options.clique, options.order, options.level)
%
[OptVal, time, stat] = solve_moment_manual(typ, obj, MomConst, LocConst, options);
end