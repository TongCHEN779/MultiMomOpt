function chordal_extension(A; turn="off")
    n = size(A, 1);
    if isequal(turn, "on")
        G = SimpleGraph(n);
        for i = 1:n-1
            for j = i+1:n
                if A[i,j] != 0
                    add_edge!(G, i, j);
                end
            end
        end
        g = gplot(G);
    end
    A = sparse(1*(A.!=0)) + (2*n + 1)*sparse(I, n, n);
    p = amd(A);
    chol = cholesky(Matrix(A[p, p]));
    clique = sparse(1*(chol.U.!=0));
    OrigIdx = sortperm(p);
    RemainIdx = [1];
    for i = 2:n
        idx = i:n;
        num = findall(clique[i, idx].!=0);
        NumOne = length(num);
        CliqueResult = sum(clique[RemainIdx, idx[num]], dims = 2);
        if isempty(findall(CliqueResult.==NumOne))
            push!(RemainIdx, i);
        end
    end
    set = clique[RemainIdx, OrigIdx];
    elem = [];
    for i = 1:size(set, 1)
        elem = vcat(elem, findall(set[i,:].!=0));
    end
    NumElem = Matrix(sum(set, dims = 2));
    cliques = Array{Any}(undef, length(NumElem));
    for i = 1:length(NumElem)
        cliques[i] = elem[sum(NumElem[1:i-1])+1:sum(NumElem[1:i])]';
    end
    G_new = SimpleGraph(n);
    for i = 1:length(NumElem)
        for j = 1:NumElem[i]-1
            for k = j+1:NumElem[i]
                add_edge!(G_new, elem[sum(NumElem[1:i-1]) + j], elem[sum(NumElem[1:i-1]) + k]);
            end
        end
    end
    A_new = adjacency_matrix(G_new);
    if isequal(turn, "on")
        g_new = gplot(G_new);
        return cliques, NumElem, A_new, g, g_new
    else isequal(turn, "off")
        return cliques, NumElem, A_new
    end
end
# vars = matread("A.mat");
# A = vars["A"];
# c,s,nA,g,ng = chordal_extension(A, turn="on");

function get_basis(n, d)
    num = binomial(n+d, d);
    basis = zeros(Int64, num, n);
    i = 0; t = 1;
    while i < d+1
        if basis[t, n] == i
            if i < d
                t = t+1;
                basis[t, 1] = i+1;
            end
            i = i+1;
        else
            j = 1;
            while basis[t, j] == 0
                j = j+1;
            end
            if j == 1
                t = t+1;
                basis[t, :] = basis[t-1, :];
                basis[t, 1] = basis[t, 1] - 1;
                basis[t, 2] = basis[t, 2] + 1;
            else
                t = t+1;
                basis[t, :] = basis[t-1, :];
                basis[t, 1] = basis[t, j] - 1;
                basis[t, j] = 0;
                basis[t, j+1] = basis[t, j+1] + 1;
            end
        end
    end
    basis = sparse(basis)
    return basis
end
# b = get_basis(2,2);

function basis2supp(basis, x)
    NumBasis = length(basis); NumVar = length(x);
    supp = Array{Int64, 2}(undef, NumBasis, NumVar);
    for i = 1:NumBasis
        for j = 1:NumVar
            supp[i,j] = MultivariatePolynomials.degree(basis[i], x[j]);
        end
    end
    supp = sparse(supp);
    return supp
end
# @polyvar x[1:3];
# basis = [x[1]; x[1]*x[3]];
# supp = basis2supp(basis, x);

function pol2supp(pol, x)
    mon = monomials(pol);
    coe = coefficients(pol);
    NumTerms = length(mon); NumVar = length(x);
    supp = Array{Float64, 2}(undef, NumTerms, 1+NumVar);
    for i = 1:NumTerms
        supp[i,1] = coe[i];
        for j = 1:NumVar
            supp[i,j+1] = MultivariatePolynomials.degree(mon[i], x[j]);
        end
    end
    supp = sparse(supp);
    return supp
end
# @polyvar x[1:3];
# pol = 1.5 + x[1] + x[2]^2*x[3];
# supp = pol2supp(pol, x);

function MomentMatrix(basis, order, SetVars, ObjQuad, model)
    if order == 1
        vars = SetVars;
        n = size(basis, 1);
        E = sparse([zeros(1+n, n-1) [2 zeros(1,n); ones(n,1) I(n)]]);
        matrix = @variable(model, [1:n+1, 1:n+1]);
        for i = 1:n+1
            if isequal(ObjQuad, "no")
                for j = i:n+1
                    new = 1;
                    for k = 1:length(vars["supp"])
                        if isequal(vars["supp"][k,:], E[j-i+1, n-i+2:2*n-i+1])
                            # @printf("%d, %d, %d\n", i,j,k)
                            matrix[j,i] = vars["var"][k];
                            matrix[i,j] = vars["var"][k];
                            new = 0;
                            break
                        end
                    end
                    if new == 1
                        vars["var"] = vcat(vars["var"], matrix[j,i]);
                        vars["supp"] = vcat(vars["supp"], E[j-i+1, n-i+2:2*n-i+1]');
                    end
                end
            elseif isequal(ObjQuad, "yes")
                # @printf("%d\n", i)
                matrix[i:n+1, i] = vars["var"][floor(Int64, (2*n+4-i)*(i-1)/2 + 1):floor(Int64, (2*n+3-i)*i/2)];
                matrix[i, i:n+1] = vars["var"][floor(Int64, (2*n+4-i)*(i-1)/2 + 1):floor(Int64, (2*n+3-i)*i/2)];
            end
        end
    elseif order >= 2
        pol = [1 zeros(1, size(basis,2))];
        vars, matrix = LocalizationMatrix(pol, basis, order, SetVars, model);
    end
    return vars, matrix
end
# n = 2; order = 2; num = binomial(n+2*order, 2*order); @polyvar x[1:num];
# basis = I(n); SetVars = x; ObjQuad = "yes";
# SetVars = Dict();
# SetVars["var"] = @variable(model, y[1:num, 1:1]);
# SetVars["supp"] = get_basis(n,2*order);
# model = Model(with_optimizer(Mosek.Optimizer));
# vars, matrix = MomentMatrix(basis, order, SetVars, ObjQuad, model);

function LocalizationMatrix(pol, basis, order, SetVars, model)
    n = size(basis, 1);
    B = get_basis(n, order); s = size(B,1); C = B*basis;
    matrix = zeros(s,s);
    de = size(pol, 1);
    vars = SetVars;
    for i = 1:s
        for j = i:s
            for k = 1:de
                new = 1;
                for z = 1:length(vars["supp"])
                    if isequal(vars["supp"][z, :], pol[k, 2:end] + C[i, :] + C[j, :])
                        E = zeros(s,s); E[i,j] = 1; E[j,i] = 1; E = sparse(E);
                        matrix = matrix + pol[k,1]*vars["var"][z]*E;
                        new = 0;
                        break
                    end
                end
                if new == 1
                    vars["var"] = vcat(vars["var"], @variable(model, [1:1, 1:1]));
                    vars["supp"] = vcat(vars["supp"], sparse(pol[k, 2:end] + C[i, :] + C[j, :])');
                    E = zeros(s,s); E[i,j] = 1; E[j,i] = 1; E = sparse(E);
                    matrix = matrix + pol[k,1]*vars["var"][length(vars["var"])]*E;
                end
            end
        end
    end
    return vars, matrix
end
# n = 2; order = 2; num = binomial(n+2*order, 2*order); @polyvar x[1:num];
# basis = I(n); SetVars = x;
# pol = [1 zeros(1, size(basis,2))];
# SetVars = Dict();
# SetVars["var"] = @variable(model, y[1:num, 1:1]);
# SetVars["supp"] = get_basis(n, 2*order);
# model = Model(with_optimizer(Mosek.Optimizer));
# vars, matrix = LocalizationMatrix(pol, basis, order, SetVars, model);
