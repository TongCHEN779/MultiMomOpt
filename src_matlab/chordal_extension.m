function [cliques, sizes] = chordal_extension(A, turn)
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
lb = size(A,1);
A = spones(A) + (2*lb+1)*speye(lb);
opts.dense=dense;
I = amd(A,opts);
R = chol(A(I,I));
Cliques = spones(R);
[~,orig_idx] = sort(I);
remainIdx = 1;
for i=2:lb
    idx = i:lb;
    one = find(Cliques(i,idx));
    noOfone = length(one);
    cliqueResult = sum(Cliques(remainIdx,idx(one)),2);
    if isempty(find(cliqueResult == noOfone,1))
        remainIdx = [remainIdx;i];
    end
end
cSet = Cliques(remainIdx,orig_idx);
[Elem,~] = find(cSet');
NoElem = full(sum(cSet,2));
cliques = cell(length(NoElem),1); sizes = [];
for i = 1:length(NoElem)
    cliques{i} = Elem(sum(NoElem(1:i-1))+1:sum(NoElem(1:i)))';
    sizes = [sizes, NoElem(i)];
end
% draw the chordal graph
if isequal(turn, 'on')
    s = []; t = [];
    for i = 1:length(NoElem)
        for j = 1:NoElem(i)-1
            for k = j+1:NoElem(i)
                s = [s; Elem(sum(NoElem(1:i-1))+j)];
                t = [t; Elem(sum(NoElem(1:i-1))+k)];
            end
        end
    end
    G = graph(s,t);
    if ismultigraph(G)
        G = simplify(G);
    end
    plot(G)
end
end