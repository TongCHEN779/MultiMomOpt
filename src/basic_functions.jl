function ChordalExtension(A; turn="off")
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

function ChordalLip(n1, n2)
    G = zeros(2*n1+n2, 2*n1+n2)
    for i = 1:n1
        for j = n1+1:n1+n2
            G[i,j] = 1; G[j,i] = 1;
        end
    end
    for i = n1+1:n1+n2
        for j = n1+n2+1:2*n1+n2
            G[i,j] = 1; G[j,i] = 1;
        end
    end
    for i = 1:n1
        for j = i:n1
            G[i,j] = 1; G[j,i] = 1;
        end
    end
    for i = n1+1:n1+n2
        G[i,i] = 1;
    end
    cliques, sizes = ChordalExtension(G)
    return cliques, sizes, G
end

function ChordalCert(n1, n2)
    # G = zeros(n1+n2, n1+n2)
    # for i = 1:n1
    #     for j = n1+1:n1+n2
    #         G[i,j] = 1; G[j,i] = 1;
    #     end
    # end
    # for i = 1:n1
    #     for j = i:n1
    #         G[i,j] = 1; G[j,i] = 1;
    #     end
    # end
    # for i = n1+1:n1+n2
    #     G[i,i] = 1
    # end
    # cliques, sizes = ChordalExtension(G)
    cliques = Array{Any}(undef, n2); sizes = Int.(zeros(n2,1));
    for i = 1:n2
        cliques[i] = vcat(collect(1:n1), n1+i)';
        sizes[i] = n1+1;
    end
    return cliques, sizes
end

function Chordal0032()
    G = zeros(100, 100)
    for i = 1:50
        for j = i:50
            G[i,j] = 1; G[j,i] = 1;
        end
    end
    for i = 1:50
        G[i,i+50] = 1; G[i+50, i] = 1
    end
    for i = 51:100
        for j = i:100
            G[i,j] = 1; G[j,i] = 1;
        end
    end
    cliques, sizes = ChordalExtension(G)
    return cliques, sizes
end

function gen_basis(n, d)
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

function heuristic(M, L, set, idx_i, depth, options; clique = [], cliques = [], siz = [], sizes = [])
    if options["clique"] == "off"
        siz = size(M, 1) - 1;
    end
    if options["method"] == "random"
        if isempty(idx_i)
            idx = sort(sample(1:siz, depth, replace=false))
        else
            idx = sort(sample(1:siz-1, depth, replace=false))
        end
    elseif options["method"] == "ordered"
        if isempty(idx_i)
            idx = collect(1:depth)
        else
            idx = collect(idx_i[1]:idx_i[1]+depth-1)
        end
    elseif options["method"] == "moment"
        norm = []
        for k = 1:siz-1
            subset = vcat(idx_i, set[k:k+options["level"]-2])
            norm = vcat(norm, sum(abs.(M[vcat(1,subset), vcat(1, subset)])))
        end
        norm_ordered = unique(sort(norm, rev=true))
        idx = []
        for k in 1:depth
            idx_orig = findall(norm.==norm_ordered[k])
            if length(idx) + length(idx_orig) <= depth
                idx = vcat(idx, idx_orig) # index of the top k-th norm
            else
                idx = vcat(idx, idx_orig[1:depth-length(idx)])
                break
            end
        end
    elseif options["method"] == "laplacian"
        norm = []
        for k = 1:siz-1
            subset = vcat(idx_i, set[k:k+options["level"]-2])
            norm = vcat(norm, sum(abs.(L[subset, subset])))
        end
        norm_ordered = unique(sort(norm, rev=true))
        idx = []
        for k = 1:depth
            idx_orig = findall(norm.==norm_ordered[k])
            if length(idx) + length(idx_orig) <= depth
                idx = vcat(idx, idx_orig) # index of the top k-th norm
            else
                idx = vcat(idx, idx_orig[1:depth-length(idx)])
                break
            end
        end
    elseif options["method"] == "mac_max"
        if options["clique"] == "off"
            error("NotImplementedError")
        end
        num = [];
        for k = 1:siz-1
            num_t = 0;
            subset = vcat(idx_i, set[k:k+options["level"]-2])
            for l = 1:length(sizes)
                if length(findall(cliques[l].-sort(subset).==0)) == options["level"]
                    num_t += 1;
                end
            end
            num = vcat(num, num_t)
        end
        num_ordered = sort(num, rev=true)
        idx = []
        for k = 1:depth
            idx_orig = findall(num.==num_ordered[k])
            if length(idx) + length(idx_orig) <= depth
                idx = vcat(idx, idx_orig)
            else
                idx = vcat(idx, idx_orig[1:depth-length(idx)])
                break
            end
        end
    elseif options["method"] == "mac_min"
        if options["clique"] == "off"
            error("NotImplementedError")
        end
        num = [];
        for k = 1:siz-1
            num_t = 0;
            subset = vcat(idx_i, set[k:k+options["level"]-2])
            for l = 1:length(sizes)
                if length(findall(cliques[l].-sort(subset).==0)) == options["level"]
                    num_t += 1;
                end
            end
            num = vcat(num, num_t)
        end
        num_ordered = sort(num)
        idx = []
        for k = 1:depth
            idx_orig = findall(num.==num_ordered[k])
            if length(idx) + length(idx_orig) <= depth
                idx = vcat(idx, idx_orig)
            else
                idx = vcat(idx, idx_orig[1:depth-length(idx)])
                break
            end
        end
    elseif options["method"] == "mac_max+moment"
        if options["clique"] == "off"
            error("NotImplementedError")
        end
        num = [];
        for k = 1:siz-1
            num_t = 0;
            subset = vcat(idx_i, set[k:k+options["level"]-2])
            for l = 1:length(sizes)
                if length(findall(cliques[l].-sort(subset).==0)) == options["level"]
                    num_t += 1;
                end
            end
            num = vcat(num, num_t)
        end
        num_ordered = sort(num, rev=true)
        idx = []
        for k = 1:depth
            idx_orig = findall(num.==num_ordered[k])
            if length(idx) + length(idx_orig) <= depth
                idx = vcat(idx, idx_orig)
            else
                norm = []
                for l = 1:length(idx_orig)
                    subset = vcat(idx_i, set[idx_orig[l]:idx_orig[l]+options["level"]-2])
                    norm = vcat(norm, sum(abs.(M[vcat(1, subset), vcat(1, subset)])))
                end
                norm_ordered = unique(sort(norm, rev=true))
                for l = 1:depth-length(idx)
                    idx_norm = findall(norm.==norm_ordered[l])
                    idx = vcat(idx, idx_orig[idx_norm])
                end
                break
            end
        end
    elseif options["method"] == "mac_min+moment"
        if options["clique"] == "off"
            error("NotImplementedError")
        end
        num = [];
        for k = 1:siz-1
            num_t = 0;
            subset = vcat(idx_i, set[k:k+options["level"]-2])
            for l = 1:length(sizes)
                if length(findall(cliques[l].-sort(subset).==0)) == options["level"]
                    num_t += 1;
                end
            end
            num = vcat(num, num_t)
        end
        num_ordered = sort(num)
        idx = []
        for k = 1:depth
            idx_orig = findall(num.==num_ordered[k])
            if length(idx) + length(idx_orig) <= depth
                idx = vcat(idx, idx_orig)
            else
                norm = []
                for l = 1:length(idx_orig)
                    subset = vcat(idx_i, set[idx_orig[l]:idx_orig[l]+options["level"]-2])
                    norm = vcat(norm, sum(abs.(M[vcat(1, subset), vcat(1, subset)])))
                end
                norm_ordered = unique(sort(norm, rev=true))
                for l = 1:depth-length(idx)
                    idx_norm = findall(norm.==norm_ordered[l])
                    idx = vcat(idx, idx_orig[idx_norm])
                end
                break
            end
        end
    elseif options["method"] == "mac_max+laplacian"
        if options["clique"] == "off"
            error("NotImplementedError")
        end
        num = [];
        for k = 1:siz-1
            num_t = 0;
            subset = vcat(idx_i, set[k:k+options["level"]-2])
            for l = 1:length(sizes)
                if length(findall(cliques[l].-sort(subset).==0)) == options["level"]
                    num_t += 1;
                end
            end
            num = vcat(num, num_t)
        end
        num_ordered = sort(num, rev=true)
        idx = []
        for k = 1:depth
            idx_orig = findall(num.==num_ordered[k])
            if length(idx) + length(idx_orig) <= depth
                idx = vcat(idx, idx_orig)
            else
                norm = []
                for l = 1:length(idx_orig)
                    subset = vcat(idx_i, set[idx_orig[l]:idx_orig[l]+options["level"]-2])
                    norm = vcat(norm, sum(abs.(L[subset, subset])))
                end
                norm_ordered = unique(sort(norm, rev=true))
                for l = 1:depth-length(idx)
                    idx_norm = findall(norm.==norm_ordered[l])
                    idx = vcat(idx, idx_orig[idx_norm])
                end
                break
            end
        end
    elseif options["method"] == "mac_min+laplacian"
        if options["clique"] == "off"
            error("NotImplementedError")
        end
        num = [];
        for k = 1:siz-1
            num_t = 0;
            subset = vcat(idx_i, set[k:k+options["level"]-2])
            for l = 1:length(sizes)
                if length(findall(cliques[l].-sort(subset).==0)) == options["level"]
                    num_t += 1;
                end
            end
            num = vcat(num, num_t)
        end
        num_ordered = sort(num)
        idx = []
        for k = 1:depth
            idx_orig = findall(num.==num_ordered[k])
            if length(idx) + length(idx_orig) <= depth
                idx = vcat(idx, idx_orig)
            else
                norm = []
                for l = 1:length(idx_orig)
                    subset = vcat(idx_i, set[idx_orig[l]:idx_orig[l]+options["level"]-2])
                    norm = vcat(norm, sum(abs.(L[subset, subset])))
                end
                norm_ordered = unique(sort(norm, rev=true))
                for l = 1:depth-length(idx)
                    idx_norm = findall(norm.==norm_ordered[l])
                    idx = vcat(idx, idx_orig[idx_norm])
                end
                break
            end
        end
    else
        error("NotImplementedError")
    end
    return idx
end

function MomentMatrix(model, basis, order, mom, vars; ObjQuad = true)
    n = size(basis, 1); # number of variables
    if order == 1 # 1st-order moment matrices
        vars = vars;
        if ObjQuad == true && isequal(options["clique"], "off") # dense 1st-order moment matrix
            matrix = @variable(model, [1:n+1, 1:n+1], Symmetric)
        elseif ObjQuad == true && isequal(options["clique"], "on") # sparse 1st-order moment matrices
            if !isassigned(vars["MomMat"], 1) # the first moment matrix
                matrix = @variable(model, [1:n+1, 1:n+1], Symmetric)
            else # check repeated monomials
                B = gen_basis(n, order)*basis; # 1st-order monomial basis
                matrix = Array{GenericAffExpr{Float64,VariableRef},2}(undef, n+1, n+1)
                # by symmetry, i and j are indices in the current moment matrix
                for i = 1:n+1
                    for j = i:n+1
                        # indices of variables: x_idxi, y_idxj; by convention, x0 = 1
                        idxi = findall(B[i,:] .== 1)
                        idxj = findall(B[j,:] .== 1)
                        for k = 1:length(vars["MomMat"])
                            if !isassigned(vars["MomMat"], k) # has searched up to the current matrix
                                break
                            elseif isempty(idxi) && isempty(idxj) # x_idxi = x_idxj = 1
                                @inbounds matrix[i,j] = 1
                                break
                            elseif isempty(idxi) && !isempty(idxj) # x_idxi = 1, x_idxj != 1
                                if sum(mom[k]["basis"][:,idxj[1]]) == 1 # check if x_idxj is in the tested clique
                                    idx = findall(mom[k]["basis"][:,idxj[1]] .== 1) # index of x_idxj in the tested clique
                                    @inbounds matrix[i,j] = vars["MomMat"][k][1, 1+idx[1]]
                                    @inbounds matrix[j,i] = matrix[i,j]
                                    break
                                else
                                    continue
                                end
                            elseif sum(mom[k]["basis"][:,idxi[1]]) == 1 && sum(mom[k]["basis"][:,idxj[1]]) == 1
                                if i == j
                                    idx = findall(mom[k]["basis"][:,idxi[1]] .== 1)
                                    @inbounds matrix[i,j] = vars["MomMat"][k][1+idx[1], 1+idx[1]]
                                elseif i != j
                                    idx1 = findall(mom[k]["basis"][:,idxi[1]] .== 1)
                                    idx2 = findall(mom[k]["basis"][:,idxj[1]] .== 1)
                                    @inbounds matrix[i,j] = vars["MomMat"][k][1+idx1[1], 1+idx2[1]]
                                    @inbounds matrix[j,i] = matrix[i,j]
                                end
                                break
                            else
                                continue
                            end
                        end
                        if !isassigned(matrix, j+(n+1)*(i-1)) # assign new variables to the unrepeated monomials
                            @inbounds matrix[i,j] = @variable(model, [1:1, 1:1])[1];
                            @inbounds matrix[j,i] = matrix[i,j];
                        end
                    end
                end
            end
        else # non-quadratic optimization
            error("NotImplementedError")
        end
    else # higher-order moment matrices
        m1 = binomial(n+1*order, 1*order); # number of monomials of degree <= order
        m2 = binomial(n+2*order, 2*order); # number of monomials of degree <= 2*order
        B = gen_basis(n, 2*order)*basis; # 2*order-degree monomial basis
        matrix = spzeros(m1, m1);
        file = matopen("Moment_and_Localizing_Matrices.mat");
        matrix_t = read(file, "moment_$(n)_$(order)"); # M = ∑y_i B_i, matrix_t is {B_i}
        close(file);
        for t = 1:m2
            idx = findall(!iszero, B[t,:]); # indices of variables
            new = 1; # new = 1, new monomial; new = 0, existed monomial
            for i = 1:length(vars["MomMat"]) # check previous moment matrices
                if !isassigned(vars["MomMat"], i)
                    break
                # check if the degree of monomial is larger than the max deg in the moment matrix
                elseif sum(B[t,:]) > 2*mom[i]["ord"]
                    continue
                elseif sum(B[t,:]) <= 2*mom[i]["ord"]
                    if isempty(idx) # constant monomial
                        new = 0;
                        @inbounds matrix = matrix .+ matrix_t[t]*1;
                        break
                    elseif sum(mom[i]["basis"][:,idx]) == length(idx) # variables appear in this clique
                        # in the 2nd degree monomial basis:
                        # the index of x_i = 1 + i, of x_i^2 = 1 + n + i (i + 1) / 2, of x_i x_j = 1 + n + i + j (j - 1) / 2
                        new = 0
                        # indices of variables in the tested clique
                        id = findall(!iszero, sum(mom[i]["basis"][:,idx], dims=2)[:])
                        if length(idx) == 1
                            if sum(B[t,:]) == 1 # x_i
                                # @printf("t=%d, i=%d, id=%d\n", t, i, id[1])
                                @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1,1+id[1]];
                                # println(matrix[1,t]==matrix_t[t][1,t]*vars["MomMat"][i][1,1+id[1]])
                            elseif sum(B[t,:]) == 2 # x_i^2
                                @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id[1],1+id[1]];
                            elseif sum(B[t,:]) == 3 # x_i^3 = x_i x_i^2
                                id_s = Int(id[1]*(id[1]+1)/2+size(mom[i]["basis"],1));
                                @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id[1],1+id_s[1]];
                            elseif sum(B[t,:]) == 4 # x_i^4 = x_i^2 x_i^2
                                id_s = Int(id[1]*(id[1]+1)/2+size(mom[i]["basis"],1));
                                @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_s[1],1+id_s[1]];
                            end
                        elseif length(idx) == 2
                            if sum(B[t,:]) == 2 # x_i x_j
                                @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id[1],1+id[2]];
                            elseif sum(B[t,:]) == 3
                                if B[t,idx[1]] == 1 # x_i x_j^2
                                    id_s = Int(id[2]*(id[2]+1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id[1],1+id_s[1]];
                                else # x_i^2 x_j
                                    id_s = Int(id[1]*(id[1]+1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_s[1],1+id[2]];
                                end
                            elseif sum(B[t,:]) == 4
                                if B[t,idx[1]] == 1 # x_i x_j^3 = x_i x_j x_j^2
                                    id_c = Int(id[1]+id[2]*(id[2]-1)/2+size(mom[i]["basis"],1));
                                    id_s = Int(id[2]*(id[2]+1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_c[1],1+id_s[1]];
                                elseif B[t,idx[1]] == 2 # x_i^2 x_j^2
                                    id_s1 = Int(id[1]*(id[1]+1)/2+size(mom[i]["basis"],1));
                                    id_s2 = Int(id[2]*(id[2]+1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_s1[1],1+id_s2[1]];
                                else # x_i^3 x_j = x_i^2 x_i x_j
                                    id_s = Int(id[1]*(id[1]+1)/2+size(mom[i]["basis"],1));
                                    id_c = Int(id[1]+id[2]*(id[2]-1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_s[1],1+id_c[1]];
                                end
                            end
                        elseif length(idx) == 3
                            if sum(B[t,:]) == 3 # x_i x_j x_k
                                id_c = Int(id[2]+id[3]*(id[3]-1)/2+size(mom[i]["basis"],1));
                                @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id[1],1+id_c[1]];
                            elseif sum(B[t,:]) == 4
                                if B[t,idx[1]] == 2 # x_i^2 x_j x_k
                                    id_s = Int(id[1]*(id[1]+1)/2+size(mom[i]["basis"],1));
                                    id_c = Int(id[2]+id[3]*(id[3]-1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_s[1],1+id_c[1]];
                                elseif B[t,idx[2]] == 2 # x_i x_j^2 x_k = x_i x_k x_j^2
                                    id_s = Int(id[2]*(id[2]+1)/2+size(mom[i]["basis"],1));
                                    id_c = Int(id[1]+id[3]*(id[3]-1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_s[1],1+id_c[1]];
                                else # x_i x_j x_k^2
                                    id_s = Int(id[3]*(id[3]+1)/2+size(mom[i]["basis"],1));
                                    id_c = Int(id[1]+id[2]*(id[2]-1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_s[1],1+id_c[1]];
                                end
                            end
                        elseif length(idx) == 4 # x_i x_j x_k x_l
                            id_c1 = Int(id[1]+id[2]*(id[2]-1)/2+size(mom[i]["basis"],1));
                            id_c2 = Int(id[3]+id[4]*(id[4]-1)/2+size(mom[i]["basis"],1));
                            @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_c1[1],1+id_c2[1]];
                        else # monomials of degree > 4
                            error("NotImplementedError")
                        end
                        break
                    else # pass to the next moment matrix
                        new = 1
                        continue
                    end
                end
            end
            if new == 1 # not existed monomial, create new variable
                var = @variable(model, [1:1, 1:1])[1];
                @inbounds matrix = matrix .+ matrix_t[t]*var;
            end
        end
    end
    matrix = Matrix(matrix);
    return matrix
end

function LocalizingMatrix(model, pol, basis, order, mom, loc, vars)
    n = size(basis, 1); # number of variables
    m1 = binomial(n+1*order,1*order); # number of monomials of degree <= order
    m2 = binomial(n+2*order,2*order); # number of monomials of degree <= 2*order
    B = gen_basis(n,2*order)*basis; # 2*order-degree monomial basis
    # M = ∑y_i B_i, matrix_t is {B_i}
    if order == 0
        matrix_t = 1;
    else
        file = matopen("Moment_and_Localizing_Matrices.mat");
        matrix_t = read(file, "moment_$(n)_$(order)");
        close(file);
    end
    de = size(pol, 1); # number of monomial multiplier
    matrix_sep = Array{Any}(undef, de); # matrices w.r.t. each monomial multiplier
    matrix_loc = 0; # localizing matrix
    for d = 1:de
        matrix = spzeros(m1,m1); # moment matrix without multiplying the monomials
        for t = 1:m2
            new = 1; # new = 1, new monomial; new = 0, existed monomial
            for i = 1:length(vars["MomMat"]) # check previous moment matrices
                C = B[t,:] .+ pol[d,2:end]
                idx = findall(!iszero, C); # indices of variables
                if !isassigned(vars["MomMat"], i)
                    break
                # check if the degree of monomial is larger than the max deg in the moment matrix
                elseif sum(C) > 2*mom[i]["ord"]
                    continue
                elseif sum(C) <= 2*mom[i]["ord"]
                    if isempty(idx) # constant monomial
                        new = 0;
                        @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1,1];
                        break
                    elseif sum(mom[i]["basis"][:,idx]) == length(idx) # variables appear in this clique
                        # in the 2nd degree monomial basis:
                        # the index of x_i = 1 + i, of x_i^2 = 1 + n + i (i + 1) / 2, of x_i x_j = 1 + n + i + j (j - 1) / 2
                        new = 0
                        # indices of variables in the tested clique
                        id = findall(!iszero, sum(mom[i]["basis"][:,idx], dims=2)[:])
                        if length(idx) == 1
                            if sum(C) == 1 # x_i
                                @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1,1+id[1]];
                            elseif sum(C) == 2 # x_i^2
                                # @printf("d=%d, t=%d, i=%d\n", d, t, i)
                                @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id[1],1+id[1]];
                            elseif sum(C) == 3 # x_i^3 = x_i x_i^2
                                # @printf("d=%d, t=%d, i=%d\n", d, t, i)
                                id_s = Int(id[1]*(id[1]+1)/2+size(mom[i]["basis"],1));
                                @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id[1],1+id_s[1]];
                            elseif sum(C) == 4 # x_i^4 = x_i^2 x_i^2
                                # @printf("d=%d, t=%d, i=%d\n", d, t, i)
                                id_s = Int(id[1]*(id[1]+1)/2+size(mom[i]["basis"],1));
                                @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_s[1],1+id_s[1]];
                            end
                        elseif length(idx) == 2
                            if sum(C) == 2 # x_i x_j
                                @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id[1],1+id[2]];
                            elseif sum(C) == 3
                                if C[idx[1]] == 1 # x_i x_j^2
                                    id_s = Int(id[2]*(id[2]+1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id[1],1+id_s[1]];
                                else # x_i^2 x_j
                                    id_s = Int(id[1]*(id[1]+1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_s[1],1+id[2]];
                                end
                            elseif sum(C) == 4
                                if C[idx[1]] == 1 # x_i x_j^3 = x_i x_j x_j^2
                                    id_c = Int(id[1]+id[2]*(id[2]-1)/2+size(mom[i]["basis"],1));
                                    id_s = Int(id[2]*(id[2]+1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_c[1],1+id_s[1]];
                                elseif C[idx[1]] == 2 # x_i^2 x_j^2
                                    id_s1 = Int(id[1]*(id[1]+1)/2+size(mom[i]["basis"],1));
                                    id_s2 = Int(id[2]*(id[2]+1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_s1[1],1+id_s2[1]];
                                else # x_i^3 x_j = x_i^2 x_i x_j
                                    id_s = Int(id[1]*(id[1]+1)/2+size(mom[i]["basis"],1));
                                    id_c = Int(id[1]+id[2]*(id[2]-1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_s[1],1+id_c[1]];
                                end
                            end
                        elseif length(idx) == 3
                            if sum(C) == 3 # x_i x_j x_k
                                id_c = Int(id[2]+id[3]*(id[3]-1)/2+size(mom[i]["basis"],1));
                                @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id[1],1+id_c[1]];
                            elseif sum(C) == 4
                                if C[idx[1]] == 2 # x_i^2 x_j x_k
                                    id_s = Int(id[1]*(id[1]+1)/2+size(mom[i]["basis"],1));
                                    id_c = Int(id[2]+id[3]*(id[3]-1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_s[1],1+id_c[1]];
                                elseif C[idx[2]] == 2 # x_i x_j^2 x_k = x_i x_k x_j^2
                                    id_s = Int(id[2]*(id[2]+1)/2+size(mom[i]["basis"],1));
                                    id_c = Int(id[1]+id[3]*(id[3]-1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_s[1],1+id_c[1]];
                                else # x_i x_j x_k^2
                                    id_s = Int(id[3]*(id[3]+1)/2+size(mom[i]["basis"],1));
                                    id_c = Int(id[1]+id[2]*(id[2]-1)/2+size(mom[i]["basis"],1));
                                    @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_s[1],1+id_c[1]];
                                end
                            end
                        elseif length(idx) == 4 # x_i x_j x_k x_l
                            id_c1 = Int(id[1]+id[2]*(id[2]-1)/2+size(mom[i]["basis"],1));
                            id_c2 = Int(id[3]+id[4]*(id[4]-1)/2+size(mom[i]["basis"],1));
                            @inbounds matrix = matrix .+ matrix_t[t]*vars["MomMat"][i][1+id_c1[1],1+id_c2[1]];
                        else
                            error("NotImplementedError")
                        end
                        break
                    else
                        new = 1
                        continue
                    end
                end
            end
            if new == 1
                for i = 1:length(vars["LocMat"]) # check previous localizing matrices
                    if !isassigned(vars["LocMat"], i)
                        break
                    else
                        for j = 1:length(vars["LocMat"][i])
                            C = B[t,:] .+ pol[d,2:end] .- loc[i]["pol"][j,2:end]
                            # 1. check if the degree of monomial is larger than the max deg in the localizing matrix
                            # 2. check if the loc monomial multiplier is in the tested monomial
                            if sum(C) > 2*loc[i]["ord"] || !isempty(findall(C .< 0))
                                continue
                            else
                                idx = findall(!iszero, C); # indices of variables
                                if isempty(idx) # constant monomial
                                    new = 0;
                                    @inbounds matrix = matrix .+ matrix_t[t]*1;
                                    break
                                elseif sum(loc[i]["basis"][:,idx]) == length(idx) # variables appear in this clique
                                    # in the 2nd degree monomial basis:
                                    # the index of x_i = 1 + i, of x_i^2 = 1 + n + i (i + 1) / 2, of x_i x_j = 1 + n + i + j (j - 1) / 2
                                    new = 0
                                    # indices of variables in the tested clique
                                    id = findall(!iszero, sum(loc[i]["basis"][:,idx], dims=2)[:])
                                    if length(idx) == 1
                                        if sum(C) == 1 # x_i
                                            @inbounds matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1,1+id[1]];
                                        elseif sum(C) == 2 # x_i^2
                                            @inbounds matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1+id[1],1+id[1]];
                                        elseif sum(C) == 3 # x_i^3 = x_i x_i^2
                                            id_s = Int(id[1]*(id[1]+1)/2+size(loc[i]["basis"],1));
                                            @inbounds matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1+id[1],1+id_s[1]];
                                        elseif sum(C) == 4 # x_i^4 = x_i^2 x_i^2
                                            id_s = Int(id[1]*(id[1]+1)/2+size(loc[i]["basis"],1));
                                            @inbounds matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1+id_s[1],1+id_s[1]];
                                        end
                                    elseif length(idx) == 2
                                        if sum(C) == 2 # x_i x_j
                                            @inbounds matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1+id[1],1+id[2]];
                                        elseif sum(C) == 3
                                            if C[idx[1]] == 1 # x_i x_j^2
                                                id_s = Int(id[2]*(id[2]+1)/2+size(loc[i]["basis"],1));
                                                matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1+id[1],1+id_s[1]];
                                            else # x_i^2 x_j
                                                id_s = Int(id[1]*(id[1]+1)/2+size(loc[i]["basis"],1));
                                                @inbounds matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1+id_s[1],1+id[2]];
                                            end
                                        elseif sum(C) == 4
                                            if C[idx[1]] == 1 # x_i x_j^3 = x_i x_j x_j^2
                                                id_c = Int(id[1]+id[2]*(id[2]-1)/2+size(loc[i]["basis"],1));
                                                id_s = Int(id[2]*(id[2]+1)/2+size(loc[i]["basis"],1));
                                                @inbounds matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1+id_c[1],1+id_s[1]];
                                            elseif C[idx[1]] == 2 # x_i^2 x_j^2
                                                id_s1 = Int(id[1]*(id[1]+1)/2+size(loc[i]["basis"],1));
                                                id_s2 = Int(id[2]*(id[2]+1)/2+size(loc[i]["basis"],1));
                                                @inbounds matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1+id_s1[1],1+id_s2[1]];
                                            else # x_i^3 x_j = x_i^2 x_i x_j
                                                id_s = Int(id[1]*(id[1]+1)/2+size(loc[i]["basis"],1));
                                                id_c = Int(id[1]+id[2]*(id[2]-1)/2+size(loc[i]["basis"],1));
                                                @inbounds matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1+id_s[1],1+id_c[1]];
                                            end
                                        end
                                    elseif length(idx) == 3
                                        if sum(C) == 3 # x_i x_j x_k
                                            id_c = Int(id[2]+id[3]*(id[3]-1)/2+size(loc[i]["basis"],1));
                                            @inbounds matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1+id[1],1+id_c[1]];
                                        elseif sum(C) == 4
                                            if C[idx[1]] == 2 # x_i^2 x_j x_k
                                                id_s = Int(id[1]*(id[1]+1)/2+size(loc[i]["basis"],1));
                                                id_c = Int(id[2]+id[3]*(id[3]-1)/2+size(loc[i]["basis"],1));
                                                @inbounds matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1+id_s[1],1+id_c[1]];
                                            elseif C[idx[2]] == 2 # x_i x_j^2 x_k = x_i x_k x_j^2
                                                id_s = Int(id[2]*(id[2]+1)/2+size(loc[i]["basis"],1));
                                                id_c = Int(id[1]+id[3]*(id[3]-1)/2+size(loc[i]["basis"],1));
                                                @inbounds matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1+id_s[1],1+id_c[1]];
                                            else # x_i x_j x_k^2
                                                id_s = Int(id[3]*(id[3]+1)/2+size(loc[i]["basis"],1));
                                                id_c = Int(id[1]+id[2]*(id[2]-1)/2+size(loc[i]["basis"],1));
                                                @inbounds matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1+id_s[1],1+id_c[1]];
                                            end
                                        end
                                    elseif length(idx) == 4 # x_i x_j x_k x_l
                                        id_c1 = Int(id[1]+id[2]*(id[2]-1)/2+size(loc[i]["basis"],1));
                                        id_c2 = Int(id[3]+id[4]*(id[4]-1)/2+size(loc[i]["basis"],1));
                                        @inbounds matrix = matrix .+ matrix_t[t]*vars["LocMat"][i][j][1+id_c1[1],1+id_c2[1]];
                                    else
                                        error("NotImplementedError")
                                    end
                                    break
                                else
                                    new = 1
                                    continue
                                end
                            end
                        end
                    end
                end
            end
            if new == 1 # not existed monomial, create new variable
                var = @variable(model, [1:1, 1:1])[1];
                @inbounds matrix = matrix .+ matrix_t[t]*var;
            end
        end
        # moment matrices w.r.t. monomials (but without multiplying monomials)
        matrix_sep[d] = matrix;
        # localizing matrices (adding up all moment matrices multiplied with monomials)
        matrix_loc = matrix_loc .+ pol[d,1]*matrix;
    end
    return matrix_loc, matrix_sep
end
