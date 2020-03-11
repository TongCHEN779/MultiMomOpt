function [cliques, sizes, newA] = chordal_extension(A, turn)
%% Chordal extension of a graph
%%
%
% INPUT:
%   A: weight matrix or adjacency matrix (double)
%   turn: if 'on' draw the chordal graph (str)
%
% OUTPUT:
%   cliques: set of maximal cliques in the chordal graph (cell)
%   sizes: sizes of each clique in the chordal graph (double)
%
%% Author: T. Chen
%%
dense = 10;
n = size(A,1);
if isequal(turn, 'on')
    subplot(1,2,1)
    G = graph;
    for i = 1:n-1
        for j = i+1:n
            if A(i,j) ~= 0
                G = addedge(G, i, j);
            end
        end
    end
    if ismultigraph(G)
        G = simplify(G);
    end
    plot(G)
    title(sprintf('Original Graph: %d vertices, %d edges', numnodes(G), numedges(G)))
    axis('square')
end
A = spones(A) + (2*n+1)*speye(n);
opts.dense=dense;
p = amd(A, opts);
ch = chol(A(p, p));
clique = spones(ch);
[~, OrigIdx] = sort(p);
RemainIdx = 1;
for i = 2:n
    idx = i:n;
    num = find(clique(i, idx));
    NumOne = length(num);
    CliqueResult = sum(clique(RemainIdx,idx(num)), 2);
    if isempty(find(CliqueResult == NumOne, 1))
        RemainIdx = [RemainIdx; i];
    end
end
set = clique(RemainIdx, OrigIdx);
[elem,~] = find(set');
NoElem = full(sum(set, 2));
cliques = cell(length(NoElem), 1); sizes = [];
for i = 1:length(NoElem)
    cliques{i} = elem(sum(NoElem(1:i-1)) + 1:sum(NoElem(1:i)))';
    sizes = [sizes, NoElem(i)];
end
% draw the chordal graph
G = graph;
for i = 1:length(NoElem)
    for j = 1:NoElem(i)-1
        for k = j+1:NoElem(i)
            G = addedge(G, elem(sum(NoElem(1:i-1))+j), elem(sum(NoElem(1:i-1))+k));
        end
    end
end
if ismultigraph(G)
    G = simplify(G);
end
newA = adjacency(G);
if isequal(turn, 'on')
    subplot(1,2,2)
    plot(G)
    title(sprintf('Chordal Graph: %d vertices, %d edges', numnodes(G), numedges(G)))
    axis('square')
end
end