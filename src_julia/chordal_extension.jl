using AMD
using LinearAlgebra
using MAT
using SparseArrays

function chordal_extension(A)
    lb = size(A, 1);
    A = sparse(1*(A.!=0)) + (2*lb + 1)*sparse(I, lb, lb);
    p = amd(A);
    R = cholesky(Matrix(A[p, p]));
    Cliques = sparse(1*(R.U.!=0));
    org_idx = sortperm(p);
    remainIdx = [1];
    for i = 2:11
        idx = i:lb;
        num = findall(Cliques[i, idx].!=0);
        noOfone = length(num);
        cliqueResult = sum(Cliques[remainIdx, idx[num]], dims = 2);
        if isempty(findall(cliqueResult.==noOfone))
            push!(remainIdx, i);
        end
    end
    cSet = Cliques[remainIdx, org_idx];
    Elem = [];
    for i = 1:size(cSet, 1)
        Elem = vcat(Elem, findall(cSet[i,:].!=0));
    end
    NoElem = Matrix(sum(cSet, dims = 2));
    cliques = Array{Any}(undef, length(NoElem));
    sizes = [];
    for i = 1:length(NoElem)
        cliques[i] = Elem[sum(NoElem[1:i-1])+1:sum(NoElem[1:i])]';
        push!(sizes, NoElem[i]);
    end
    return cliques, sizes
end

vars = matread("A.mat");
A = vars["A"];
c, s = chordal_extension(A);
