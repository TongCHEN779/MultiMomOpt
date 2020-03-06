function [OptVal, time, stat] = solve_moment_maxcut(A, options, order, level)
%% Solve MAX-CUT problem using rounding technique
%%
%
% INPUT:
%   A: adjacency matrix (double)
%   options: solver settings (cf. sdpsettings for YALMIP)
%   order: relaxation order
%   level: level of the subhierarchy
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
n = size(A, 1); L = sparse(diag(A*ones(n,1))) - A; typ = 'max';
%
obj = []; E = sparse([eye(n,n); eye(n,n)]);
for i = 1:n
    for j = i:n
        if i == j
            obj = [obj; [L(i,j)/4, E(i,:)+E(j,:)]];
        elseif i < j
            obj = [obj; [L(i,j)/2, E(i,:)+E(j,:)]];
        end
    end
end
%
if nargin == 2 || level == 0 % Shor's relaxation
    MomConst{1}.basis = sparse(eye(n,n)); MomConst{1}.ord = 1;
    %
    for i = 1:n
        LocConst{i}.pol = sparse([1, zeros(1,n); -1, 2*E(i,:)]);
        LocConst{i}.basis = sparse(eye(n,n));
        LocConst{i}.typ = '==';
        LocConst{i}.ord = 0;
    end
elseif nargin > 2 && level > 0
    [cliques, sizes] = chordal_extension(A, 'off');
    %
    MomConst{1}.basis = sparse(eye(n,n)); MomConst{1}.ord = 1;
    NumMom = 1;
    for i = 1:length(cliques)
        clique = [cliques{i}, cliques{i}];
        if level < sizes(i)
            for j = 1:sizes(i)
                NumMom = NumMom + 1;
                MomConst{NumMom}.basis = sparse(E(clique(j:j+level-1), :));
                MomConst{NumMom}.ord = order;
            end
        else
            NumMom = NumMom + 1;
            MomConst{NumMom}.basis = sparse(E(1:size(i), :));
            MomConst{NumMom}.ord = order;
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
        if level < sizes(idx)
            for j = 1:sizes(idx)
                NumLoc = NumLoc + 1;
                LocConst{NumLoc}.pol = sparse([1, zeros(1,n); -1, 2*E(i,:)]);
                LocConst{NumLoc}.basis = sparse(E(clique(j:j+level-1), :));
                LocConst{NumLoc}.typ = '==';
                LocConst{NumLoc}.ord = order - 1;
            end
        else
            NumLoc = NumLoc + 1;
            LocConst{NumLoc}.pol = sparse([1, zeros(1,n); -1, 2*E(i,:)]);
            LocConst{NumLoc}.basis = sparse(E(1:sizes(idx), :));
            LocConst{NumLoc}.typ = '==';
            LocConst{NumLoc}.ord = order - 1;
        end
    end
end
%
[OptVal, time, stat] = solve_moment_manual(typ, obj, MomConst, LocConst, options);
end