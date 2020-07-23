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
# vars = matread("A.mat");
# A = vars["A"];
# c,s,nA,g,ng = ChordalExtension(A, turn="on");

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
# b = gen_basis(2,2);

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

function MomentMatrix(model, basis, order, SetVars; ObjQuad = true)
    n = size(basis, 1);
    # println(n)
    if order == 1
        vars = SetVars;
        E = sparse([zeros(1+n, n-1) [2 zeros(1,n); ones(n,1) I(n)]]);
        matrix = Array{GenericAffExpr{Float64,VariableRef},2}(undef, n+1, n+1); #@variable(model, [1:n+1, 1:n+1]);
        # for i = 1:n+1
        #     if ObjQuad == false
        #         for j = i:n+1
        #             new = 1;
        #             for k = 1:length(vars["var"])
        #                 if isequal(vars["supp"][k,:], E[j-i+1, n-i+2:2*n-i+1])
        #                     # @printf("%d, %d, %d\n", i,j,k)
        #                     matrix[j,i] = vars["var"][k];
        #                     matrix[i,j] = vars["var"][k];
        #                     new = 0;
        #                     break
        #                 end
        #             end
        #             if new == 1
        #                 vars["var"] = vcat(vars["var"], matrix[j,i]);
        #                 vars["supp"] = vcat(vars["supp"], E[j-i+1, n-i+2:2*n-i+1]');
        #             end
        #         end
        #     elseif ObjQuad == true # objective support must be set in matrix order (not lexicographical order)
        #         # @printf("%d\n", i)
        #         matrix[i:n+1, i] = vars["var"][floor(Int64, (2*n+4-i)*(i-1)/2 + 1):floor(Int64, (2*n+3-i)*i/2)];
        #         matrix[i, i:n+1] = vars["var"][floor(Int64, (2*n+4-i)*(i-1)/2 + 1):floor(Int64, (2*n+3-i)*i/2)];
        #         # matrix[i:n+1, i] = @variable(model, [1:n+2-i]);
        #         # matrix[i, i:n+1] = matrix[i:n+1, i];
        #     end
        # end
        if ObjQuad == true
            matrix = @variable(model, [1:n+1, 1:n+1], Symmetric)
        end
    elseif order == 2 && n <= 20
        m1 = binomial(n+2*1, 2*1); m2 = binomial(n+2*2, 2*2);
        B = gen_basis(n, 2*2)*basis;
        file = matopen("Moment_and_Localization_Matrices.mat");
        matrix_temp = read(file, "moment_$(n)_2");
        close(file);
        matrix = spzeros(m1, m1);
        SetVars["var"] = Array{GenericAffExpr{Float64,VariableRef},2}(undef, m2-m1, 1)
        SetVars["supp"] = spzeros(m2-m1, size(basis, 2))
        idx = 0;
        for t = 1:m1
            idx = findall(!iszero, B[t,:]);
            if length(idx) == 0
                matrix = matrix + matrix_temp[t]*SetVars["MomMat"][1][1,1];
            elseif length(idx) == 1 && B[t,idx[1]] == 1
                matrix = matrix + matrix_temp[t]*SetVars["MomMat"][1][1,1+idx[1]];
            elseif length(idx) == 1 && B[t,idx[1]] == 2
                matrix = matrix + matrix_temp[t]*SetVars["MomMat"][1][1+idx[1],1+idx[1]];
            elseif length(idx) == 2
                matrix = matrix + matrix_temp[t]*SetVars["MomMat"][1][1+idx[1],1+idx[2]];
            end
        end
        # K = collect(1:n)
        for i = 2:length(SetVars["MomMat"])
            if isempty(SetVars["MomMat"][i])
                break
            else
                n_orig = size(SetVars["MomConst"][i], 1)
                m1_orig = binomial(n_orig+2*1, 2*1); m2_orig = binomial(n_orig+2*2, 2*2);
                set_orig = []; set_rep = []; #set_nonrep = [];
                for t = 1:n
                    rep = true
                    for k = 1:n_orig
                        if isequal(SetVars["MomConst"][i][k,:]', basis[t,:]')
                            set_orig = vcat(set_orig, k); set_rep = vcat(set_rep, t); #global set_nonrep = vcat(set_nonrep, 0);
                            # rep = true
                            break
                        # else
                        #     rep = false
                        end
                    end
                    # if !rep
                    #     # global set_orig = vcat(set_orig, 0); global set_rep = vcat(set_rep, 0);
                    #     global set_nonrep = vcat(set_nonrep, t);
                    # end
                end
                if isempty(set_rep)
                    continue
                end
                basis_orig = copy(SetVars["MomConst"][i]); basis_orig = Int.(basis_orig) #basis_nonrep = copy(basis);
                basis_rep = copy(basis); basis_rep = Int.(basis_rep)
                # compt = 0; B_index = []
                for j = 1:n_orig
                    if isempty(findall(iszero,set_orig.-j))
                        # global basis_orig = basis_orig[1:end .!= j-compt, 1:end]
                        # global compt = compt + 1
                        basis_orig[j,:] .= 4;
                    # else
                    #     global basis_nonrep = Matrix(basis_nonrep)
                    #     basis_nonrep[j,:] .= 0;
                    #     basis_nonrep = sparse(basis_nonrep)
                    end
                end
                for j = 1:n
                    if isempty(findall(iszero,set_rep.-j))
                        basis_rep[j,:] .= 4;
                    end
                end
                # println(size(gen_basis(n, 2*2)));println(size(basis_orig))
                B_orig = gen_basis(n_orig, 2*2)*basis_orig;
                # for j = 1:m2
                #     if sum(B_orig[j,:]) > 4
                #         B_orig[j,:] .= 0
                #     end
                # end
                B_orig = B_orig[m1_orig+1:m2_orig, 1:end]; u = 1; compt = 0;
                B_index = []
                while u <= m2_orig-m1_orig
                    if sum(B_orig[u-compt,:]) > 4
                        B_orig = B_orig[1:end .!= u-compt, 1:end]
                        compt = compt + 1;
                    else
                        B_index = vcat(B_index, u)
                    end
                    u = u+1;
                end
                B_rep = gen_basis(n, 2*2)*basis_rep;
                # B_rep = B_rep[m1+1:m2, 1:end]
                for j = 1:m2
                    if sum(B_rep[j,:]) > 4
                        B_rep[j,:] .= 0
                    end
                end
                B = B .- B_rep; B = B .* (B .>= 0);
                # B_nonrep = sparse(gen_basis(n, 2*1)*Matrix(basis)) .- B_rep;
                # B_nonrep = sparse(gen_basis(n, 2*2)*Matrix(basis_nonrep));
                compt = 0;
                for t = m1+1:m2
                    if sum(B_rep[t,:]) !== 0 && !iszero(matrix[findall(!iszero, matrix_temp[t])[1]])
                        compt = compt + 1;
                    elseif sum(B_rep[t,:]) !== 0 && iszero(matrix[findall(!iszero, matrix_temp[t])[1]])
                        # println(t)
                        compt = compt + 1;
                        SetVars["var"][t-m1] = SetVars["MomVar"][i][B_index[compt]];
                        # SetVars["var"][t-m1] = SetVars["MomVar"][i][t-m1];
                        SetVars["supp"][t-m1,:] .= B_rep[t,:];
                        matrix = matrix + matrix_temp[t]*SetVars["var"][t-m1];
                    end
                end
                # filter!(e->e∉set_rep, K)
                # if !isempty(SetVars["MomMat"][11]) && isempty(SetVars["MomMat"][12]) && i == 2
                #     println(Matrix(matrix))
                # end
            end
        end
        if sum(B) != 0
            for t = m1+1:m2
                if sum(B[t,:]) !== 0 #&& iszero(matrix[findall(!iszero, matrix_temp[t])[1]])
                    # println(t)
                    SetVars["var"][t-m1] = @variable(model, [1:1, 1:1])[1];
                    SetVars["supp"][t-m1,:] .= B[t,:];
                    matrix = matrix + matrix_temp[t]*SetVars["var"][t-m1]
                end
                # @printf("%d\n", t)
                # idx = 0; new = 1; p = size(basis, 2)
                # for k = 1:length(SetVars["var"])
                #     if isequal(SetVars["supp"][k,:]', B[t,:]')
                #         # @printf("%d %d\n", t,k)
                #         idx = k; new = 0;
                #         break
                #     end
                # end
                # if new == 0
                #     matrix = matrix + matrix_temp[t]*SetVars["var"][idx];
                # elseif new == 1
                #     if isempty(SetVars["var"])
                #         SetVars["var"] = @variable(model, [1:1, 1:1]);
                #         SetVars["supp"] = B[t,:]';
                #     else
                #         SetVars["var"] = vcat(SetVars["var"], @variable(model, [1:1, 1:1]));
                #         SetVars["supp"] = vcat(SetVars["supp"], B[t,:]');
                #     end
                #     matrix = matrix + matrix_temp[t]*SetVars["var"][length(SetVars["var"])]
                # end
            end
        end
        matrix = Matrix(matrix); #MomMat = @variable(model, [1:m1, 1:m1], Symmetric)
        # @constraint(model, MomMat.==matrix)
        # if isequal(ObjQuad, true)
        #     vars = SetVars;
        #     for t = 1:m1
        #         # idx = 0;
        #         # for k = 1:length(vars["var"])
        #         #     if isequal(vars["supp"][k,:]', B[t,:]')
        #         #         # @printf("%d %d\n", t,k)
        #         #         idx = k;
        #         #         matrix = matrix + matrix_temp[t]*vars["var"][idx];
        #         #         break
        #         #     end
        #         # end
        #         idx = findall(!iszero, Matrix(B)[t,:]);
        #         if length(idx) == 0
        #             matrix = matrix + matrix_temp[t]*vars["Mat"][1,1];
        #         elseif length(idx) == 1
        #             matrix = matrix + matrix_temp[t]*vars["Mat"][1+idx[1],1+idx[1]];
        #         elseif length(idx) == 2
        #             matrix = matrix + matrix_temp[t]*vars["Mat"][1+idx[1],1+idx[2]];
        #         end
        #     end
        #     for t = m1+1:m2
        #         # @printf("%d\n", t)
        #         vars["var"] = vcat(vars["var"], @variable(model, [1:1, 1:1]));
        #         vars["supp"] = vcat(vars["supp"], B[t,:]');
        #         matrix = matrix + matrix_temp[t]*vars["var"][length(vars["var"])]
        #     end
        #     matrix = Matrix(matrix);
        # elseif isequal(ObjQuad, false)
        #     vars = SetVars;
        #     for t = 1:m2
        #         idx = 0; new = 1;
        #         for k = 1:length(vars["var"])
        #             if isequal(vars["supp"][k,:]', B[t,:]')
        #                 # @printf("%d %d\n", t,k)
        #                 idx = k; new = 0;
        #                 break
        #             end
        #         end
        #         if new == 0
        #             matrix = matrix + matrix_temp[t]*vars["var"][idx];
        #         elseif new == 1
        #             vars["var"] = vcat(vars["var"], @variable(model, [1:1, 1:1]));
        #             vars["supp"] = vcat(vars["supp"], B[t,:]');
        #             matrix = matrix + matrix_temp[t]*vars["var"][length(vars["var"])]
        #         end
        #     end
        #     matrix = Matrix(matrix);
        # end

        # vars = SetVars;
        # for t = 1:m2
        #     idx = 0; new = 1;
        #     for k = 1:length(vars["var"])
        #         if isequal(vars["supp"][k,:]', B[t,:]')
        #             # @printf("%d %d\n", t,k)
        #             idx = k; new = 0;
        #             break
        #         end
        #     end
        #     if new == 0
        #         matrix = matrix + matrix_temp[t]*vars["var"][idx];
        #     elseif new == 1
        #         vars["var"] = vcat(vars["var"], @variable(model, [1:1, 1:1]));
        #         vars["supp"] = vcat(vars["supp"], B[t,:]');
        #         matrix = matrix + matrix_temp[t]*vars["var"][t]
        #     end
        # end
        # matrix = Matrix(matrix);
    elseif (order == 2 && n > 20) || order >= 3
        pol = [1 zeros(1, size(basis,2))];
        matrix = LocalizationMatrix(model, pol, basis, order, SetVars);
    end
    return matrix
end
# start = time();
# n = 1; order = 2; num = binomial(20+2, 2);
# basis = [1 zeros(1,19)];
# SetVars = Dict();
# model = Model(with_optimizer(Mosek.Optimizer));
# SetVars["var"] = @variable(model, [1:num, 1:1]);;
# SetVars["supp"] = gen_basis(20,2);
# vars, matrix = MomentMatrix(model, basis, order, SetVars, ObjQuad = "no");
# elapsed = time() - start;

function LocalizationMatrix(model, pol, basis, order, SetVars)
    n = size(basis, 1);
    m1 = binomial(n+order,order); m2 = binomial(n+order*2,order*2); m3 = binomial(n+order*4,order*4);
    idx_rep = []
    for t = 1:n
        idx_rep = vcat(idx_rep, findall(!iszero, basis[t,:]))
    end
    matrix = spzeros(m1,m1);
    de = size(pol, 1);
    SetVars["TempLocMat"] = Array{Any}(undef, de);
    for i = 1:de
        SetVars["TempLocMat"][i] = []
    end
    if order == 0
        for k = 1:de
            idx = findall(!iszero, pol[k,2:end])
            if isempty(idx)
                matrix = matrix .+ pol[k,1]*SetVars["MomMat"][1][1,1]
            elseif sum(pol[k,2:end]) == 1
                matrix = matrix .+ pol[k,1]*SetVars["MomMat"][1][1+idx[1],1]
            elseif sum(pol[k,2:end]) == 2 && length(idx) == 2
                matrix = matrix .+ pol[k,1]*SetVars["MomMat"][1][1+idx[1],1+idx[2]]
            elseif sum(pol[k,2:end]) == 2 && length(idx) == 1
                matrix = matrix .+ pol[k,1]*SetVars["MomMat"][1][1+idx[1],1+idx[1]]
            end
        end
    else
        file = matopen("Moment_and_Localization_Matrices.mat");
        matrix_temp = read(file, "moment_$(n)_$(order)");
        close(file);
        for k = 1:de
            B = gen_basis(n, order)*basis; C = gen_basis(n, order*2)*basis; D = gen_basis(n, order*4)*basis;
            # print("iter: "); println(k)
            # k = k+1
            new = 1;
            for i = 1:size(C,1)
                if sum(abs.(pol[k, 2:end] .- C[i,:])) == 0
                    # @printf("i: %d\n",i)
                    new = 0
                    break
                end
            end
            # new
            if new == 0
                if sum(pol[k, 2:end]) == 0
                    idx = 1
                    for i = 1:size(basis,1)
                        idx = vcat(idx, findall(!iszero, basis[i,:])[:].+1)
                    end
                    matrix_t = SetVars["MomMat"][1][idx[:], idx[:]]
                else
                    matrix_t = spzeros(m1,m1);
                    for u = 1:size(C,1)
                        for v = 1:size(D,1)
                            if isequal(C[u,:]+pol[k,2:end], D[v,:]) #&& iszero(matrix_t[findall(!iszero, matrix_temp[u])[1]])
                                idx = findall(!iszero, D[v,:])
                                # println(v); println(D[v,:])
                                if sum(D[v,:]) == 1
                                    matrix_t = matrix_t + matrix_temp[u]*SetVars["MomMat"][1][1+idx[1],1];
                                elseif sum(D[v,:]) == 2 && length(idx) == 2
                                    matrix_t = matrix_t + matrix_temp[u]*SetVars["MomMat"][1][1+idx[1],1+idx[2]];
                                elseif sum(D[v,:]) == 2 && length(idx) == 1
                                    matrix_t = matrix_t + matrix_temp[u]*SetVars["MomMat"][1][1+idx[1],1+idx[1]];
                                else
                                    for j = 2:length(SetVars["MomMat"])
                                        if isequal(size(SetVars["MomConst"][j],1), size(basis,1)) && sum(abs.(SetVars["MomConst"][j] .- basis)) == 0
                                            # println([j v size(C,1)]);
                                            # println(SetVars["MomVar"][s])
                                            matrix_t = matrix_t + matrix_temp[u]*SetVars["MomVar"][j][v-size(C,1)];
                                            # @printf("j: %d\n", j)
                                            # s = j
                                            # @printf("s: %d\n", s)
                                            break
                                        end
                                    end
                                end
                                break
                            end
                        end
                    end
                end
                SetVars["TempLocMat"][k] = Matrix(matrix_t)
                matrix = matrix + pol[k,1]*matrix_t
                # println(matrix[1,1])
                # matrix[1,1]
            else
                matrix_t = spzeros(m1,m1);
                for i = 2:length(SetVars["MomMat"])
                    # i = i+1
                    idx = findall(!iszero, pol[k,2:end])
                    idx_orig = []
                    for t = 1:size(SetVars["MomConst"][i], 1)
                        idx_orig = vcat(idx_orig, findall(!iszero, SetVars["MomConst"][i][t,:]))
                    end
                    mem = 0
                    for u in idx
                        for v in idx_orig
                            if u == v
                                mem = 1
                                # println(mem)
                                break
                            else
                                mem = 0
                            end
                        end
                        if mem == 0
                            break
                        end
                        # println([k i mem])
                    end
                    # if k == 5 && i == 7
                    #     println(idx);println(idx_orig);println(mem)
                    # end
                    # mem
                    # if mem==1
                    #     println([i mem])
                    # end
                    # @printf("mem: %d\n", mem)
                    if mem == 1
                        C_t = gen_basis(size(SetVars["MomConst"][i],1), order*4)*SetVars["MomConst"][i];
                        # println(Matrix(C_t));println(Matrix(C));println(Matrix(pol[k,2:end]))
                        # println([m1 m2])
                        mem = 0
                        for t = 1:m2
                            # t = t+1
                            if (t == 1 || (t != 1 && sum(C[t,:]) != 0)) && iszero(matrix_t[findall(!iszero, matrix_temp[t])[1]])
                                for s = 1:binomial(size(SetVars["MomConst"][i],1)+order*2,order*2)
                                    # s = s+1
                                    if isequal(C_t[s,:], pol[k,2:end].+C[t,:])
                                        mem = 1
                                        idx = findall(!iszero, C_t[s,:])
                                        # println(s); println(C_t[s,:])
                                        # println(idx)
                                        for u = 1:length(idx)
                                            for v = 2:m1
                                                if C_t[v,idx[u]] == 1
                                                    idx[u] = v
                                                    break
                                                end
                                            end
                                        end
                                        # println(idx)
                                        if sum(C_t[s,:]) == 1
                                            matrix_t = matrix_t + matrix_temp[t]*SetVars["MomMat"][i][idx[1], 1];
                                        elseif sum(C_t[s,:]) == 2 && length(idx) == 2
                                            matrix_t = matrix_t + matrix_temp[t]*SetVars["MomMat"][i][idx[1], idx[2]];
                                        elseif sum(C_t[s,:]) == 2 && length(idx) == 1
                                            matrix_t = matrix_t + matrix_temp[t]*SetVars["MomMat"][i][idx[1], idx[1]];
                                        end
                                        break
                                    else
                                        mem = 0
                                    end
                                end
                                if mem == 0
                                    for s = binomial(size(SetVars["MomConst"][i],1)+order*2,order*2)+1:size(C_t,1)
                                        if isequal(C_t[s,:], pol[k,2:end].+C[t,:])
                                            mem = 1
                                            # println(s); println(C_t[s,:])
                                            matrix_t = matrix_t + matrix_temp[t]*SetVars["MomVar"][i][s-binomial(size(SetVars["MomConst"][i],1)+order*2,order*2)]
                                            # println(mem)
                                            break
                                        end
                                    end
                                end
                                if mem == 1
                                    C[t,:] .= 0
                                end
                            end
                        end
                    end
                    # idx = findall(!iszero, pol[k,2:end])
                    # s1 = size(MomConst[i]["basis"],1); s2 = size(MomConst[i]["basis"],2);
                    # mem = 1; mul = MomConst[i]["basis"]
                    # for j = 1:length(idx)
                    #     mul = abs.(MomConst[i]["basis"] .- [zeros(s1,idx[j]-1) ones(s1,1) zeros(s1,s2-idx[j])]) * ones(s2,1)
                    #     if isempty(findall(iszero, mul))
                    #         mem = 0
                    #         break
                    #     else
                    #         mem = 1
                    #     end
                    # end
                    # # mem
                    # if mem == 1
                    #     set_orig = []; set_rep = [];
                    #     for j = 1:size(MomConst[i]["basis"], 1)
                    #         for t = 1:n
                    #             if isequal(MomConst[i]["basis"][j,:]', basis[t,:]')
                    #                 set_orig = vcat(set_orig, j); set_rep = vcat(set_rep, t);
                    #                 break
                    #             end
                    #         end
                    #         # if isequal(SetVars["MomConst"][i][j,:]', pol[k,2:end]')
                    #         #     global set_orig = vcat(set_orig, j);
                    #         # end
                    #     end
                    #     basis_orig = copy(MomConst[i]["basis"]); basis_orig = Int.(basis_orig)
                    #     basis_rep = copy(basis); basis_rep = Int.(basis_rep)
                    #     for j = 1:size(gen_basis(n, order*2)*basis_orig,1)
                    #         if isequal((gen_basis(n, order*2)*basis_orig)[j,:], pol[k,2:end])
                    #             idx_pol = j
                    #             break
                    #         end
                    #     end
                    #     # idx_pol
                    #     if isempty(set_rep)
                    #         for j = 1:n
                    #             basis_orig[j,:] .= 4
                    #             basis_rep[j,:] .= 4
                    #         end
                    #     else
                    #         for j = 1:n
                    #             if isempty(findall(iszero,set_orig.-j))
                    #                 basis_orig[j,:] .= 4;
                    #             end
                    #             if isempty(findall(iszero,set_rep.-j))
                    #                 basis_rep[j,:] .= 4;
                    #             end
                    #         end
                    #     end
                    #     C_orig = sparse(gen_basis(n, order*2)*Matrix(basis_orig));
                    #     # C_orig = C_orig[m2+1:m3, 1:end];
                    #     u = 1; compt = 0;
                    #     C_index = []
                    #     while u <= m2
                    #         if sum(C_orig[u-compt,:]) > 4
                    #             C_orig = C_orig[1:end .!= u-compt, 1:end]
                    #             compt = compt + 1;
                    #         else
                    #             C_index = vcat(C_index, u)
                    #         end
                    #         u = u+1;
                    #     end
                    #     C_rep = sparse(gen_basis(n, order*2)*Matrix(basis_rep));
                    #     for j = 1:m2
                    #         if sum(C_rep[j,:]) > 4
                    #             C_rep[j,:] .= 0
                    #         end
                    #     end
                    #     C = C .- C_rep; C = C .* (C .>= 0);
                    #     compt = 1;
                    #     # println([i,idx_pol, C_index[1]])
                    #     if iszero(matrix_t[findall(!iszero, matrix_temp[1])[1]])
                    #         matrix_t = matrix_t + matrix_temp[1]*SetVars["MomMat"][i][idx_pol, C_index[1]];
                    #     end
                    #     for t = 2:m2
                    #         if sum(C_rep[t,:]) != 0 && iszero(matrix_t[findall(!iszero, matrix_temp[t])[1]])
                    #             # println(t)
                    #             compt = compt + 1;
                    #             # SetVars["var"][t-m1] = SetVars["MomVar"][i][C_index[compt]];
                    #             # SetVars["supp"][t-m1,:] .= B_rep[t,:];
                    #             matrix_t = matrix_t + matrix_temp[t]*SetVars["MomMat"][i][idx_pol, C_index[compt]];
                    #         end
                    #     end
                    # end
                end
                # if k == 10
                #     print("Mom: "); println(isequal(matrix_t[1,1], SetVars["TempLocMat"][5][1,1]))
                # elseif k == 11
                #     print("Mom: "); println(isequal(matrix_t[1,1], SetVars["TempLocMat"][6][1,1]))
                # elseif k == 12
                #     print("Mom: "); println(isequal(matrix_t[1,1], SetVars["TempLocMat"][7][1,1]))
                # end
                # print("Mom: ");println(matrix_t[1,1])
                for i = 1:length(SetVars["LocConst"])
                    if isempty(SetVars["LocMat"][i][1])
                        break
                    end
                    if sum(C) != 0
                        idx_orig = []
                        for t = 1:size(SetVars["LocConst"][i], 1)
                            idx_orig = vcat(idx_orig, findall(!iszero, SetVars["LocConst"][i][t,:]))
                        end
                        for j = 1:size(SetVars["LocPol"][i],1)
                            if isempty(SetVars["LocMat"][i][j])
                                break
                            else
                                idx = findall(!iszero, pol[k,2:end])
                                idx_t = findall(!iszero, SetVars["LocPol"][i][j,2:end])
                                mem = 0
                                for u in idx
                                    for v in idx_t
                                        if u == v
                                            mem = 1
                                            break
                                        else
                                            mem = 0
                                        end
                                    end
                                    if mem == 0
                                        for v in idx_orig
                                            if u == v
                                                mem = 1
                                                break
                                            else
                                                mem = 0
                                            end
                                        end
                                    end
                                    if mem == 0
                                        break
                                    end
                                end
                                if mem == 1
                                    mem = 0
                                    for u in idx_t
                                        for v in idx
                                            if u == v
                                                mem = 1
                                                break
                                            else
                                                mem = 0
                                            end
                                        end
                                        if mem == 0
                                            for v in idx_rep
                                                if u == v
                                                    mem = 1
                                                    break
                                                else
                                                    mem = 0
                                                end
                                            end
                                        end
                                        if mem == 0
                                            break
                                        end
                                    end
                                end
                                if mem == 1
                                    mem = 0
                                    C_t = gen_basis(size(SetVars["LocConst"][i],1), order*2)*SetVars["LocConst"][i];
                                    # println(Matrix(C_t))
                                    # println([m1 m2])
                                    for t = 1:m2
                                        if (t == 1 || (t != 1 && sum(C[t,:]) != 0)) && iszero(matrix_t[findall(!iszero, matrix_temp[t])[1]])
                                            for s = 1:size(C_t,1)
                                                if isequal(SetVars["LocPol"][i][j,2:end].+C_t[s,:], pol[k,2:end].+C[t,:])
                                                    mem = 1
                                                    idx = findall(!iszero, C_t[s,:])
                                                    # println(idx)
                                                    for u = 1:length(idx)
                                                        for v = 2:m1
                                                            if C_t[v,idx[u]] == 1
                                                                idx[u] = v
                                                                break
                                                            end
                                                        end
                                                    end
                                                    # println(idx)
                                                    if sum(C_t[s,:]) == 1
                                                        matrix_t = matrix_t + matrix_temp[t]*SetVars["LocMat"][i][j][idx[1], 1];
                                                    elseif sum(C_t[s,:]) == 2 && length(idx) == 2
                                                        matrix_t = matrix_t + matrix_temp[t]*SetVars["LocMat"][i][j][idx[1], idx[2]];
                                                    elseif sum(C_t[s,:]) == 2 && length(idx) == 1
                                                        matrix_t = matrix_t + matrix_temp[t]*SetVars["LocMat"][i][j][idx[1], idx[1]];
                                                    end
                                                    break
                                                else
                                                    mem = 0
                                                end
                                            end
                                            if mem == 0
                                                C[t,:] .= 0
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                # print("Loc: ");println(matrix_t[1,1])
                if !isempty(SetVars["TempLocMat"][1])
                    if sum(C) != 0
                        for j = 1:size(pol,1)
                            if isempty(SetVars["TempLocMat"][j])
                                break
                            else
                                idx = findall(!iszero, pol[k,2:end])
                                idx_t = findall(!iszero, pol[j,2:end])
                                mem = 0
                                for u in idx
                                    for v in idx_t
                                        if u == v
                                            mem = 1
                                            break
                                        else
                                            mem = 0
                                        end
                                    end
                                    if mem == 0
                                        for v in idx_rep
                                            if u == v
                                                mem = 1
                                                break
                                            else
                                                mem = 0
                                            end
                                        end
                                    end
                                    if mem == 0
                                        break
                                    end
                                end
                                if mem == 1
                                    mem = 0
                                    for u in idx_t
                                        for v in idx
                                            if u == v
                                                mem = 1
                                                break
                                            else
                                                mem = 0
                                            end
                                        end
                                        if mem == 0
                                            for v in idx_rep
                                                if u == v
                                                    mem = 1
                                                    break
                                                else
                                                    mem = 0
                                                end
                                            end
                                        end
                                        if mem == 0
                                            break
                                        end
                                    end
                                end
                                if mem == 1
                                    # println(j)
                                    mem = 0
                                    C_t = gen_basis(size(basis,1), order*2)*basis;
                                    for t = 1:m2
                                        if (t == 1 || (t != 1 && sum(C[t,:]) != 0)) && iszero(matrix_t[findall(!iszero, matrix_temp[t])[1]])
                                            for s = 1:size(C_t,1)
                                                if isequal(pol[j,2:end].+C_t[s,:], pol[k,2:end].+C[t,:])
                                                    mem = 1
                                                    idx = findall(!iszero, C_t[s,:])
                                                    for u = 1:length(idx)
                                                        for v = 2:m1
                                                            if C_t[v,idx[u]] == 1
                                                                idx[u] = v
                                                                break
                                                            end
                                                        end
                                                    end
                                                    if sum(C_t[s,:]) == 1
                                                        matrix_t = matrix_t + matrix_temp[t]*SetVars["TempLocMat"][j][idx[1], 1];
                                                    elseif sum(C_t[s,:]) == 2 && length(idx) == 2
                                                        matrix_t = matrix_t + matrix_temp[t]*SetVars["TempLocMat"][j][idx[1], idx[2]];
                                                    elseif sum(C_t[s,:]) == 2 && length(idx) == 1
                                                        matrix_t = matrix_t + matrix_temp[t]*SetVars["TempLocMat"][j][idx[1], idx[1]];
                                                    end
                                                    break
                                                else
                                                    mem = 0
                                                end
                                            end
                                            if mem == 0
                                                C[t,:] .= 0
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                # print("Tem: ");println(matrix_t[1,1])
                if sum(C) != 0
                    for t = 1:m2
                        if (t == 1 || (t != 1 && sum(C[t,:]) != 0)) && iszero(matrix_t[findall(!iszero, matrix_temp[t])[1]])
                            matrix_t = matrix_t + matrix_temp[t]*@variable(model, [1:1, 1:1])[1];
                        end
                    end
                end
                # print("Dif: ");println(matrix_t[1,1])
                # print(k);print(" ");println(pol[k,1]*matrix_t[1,1])
                SetVars["TempLocMat"][k] = Matrix(matrix_t)
                matrix = matrix + pol[k,1]*matrix_t;
                # println(matrix[1,1])
                # matrix[1,1]
            end
        end
    end
        # for i = 1:s
        #     for j = i:s
        #         for k = 1:de
        #             new = 1; p = size(basis, 2); v = pol[k, 2:end] + B[i, :] + B[j, :]
        #             idx = findall(!iszero, v);
        #             E = zeros(s,s); E[i,j] = 1; E[j,i] = 1; E = sparse(E);
        #             if length(idx) == 0
        #                 matrix = matrix + pol[k,1]*SetVars["MomMat"][1][1,1]*E;
        #             elseif length(idx) == 1 && v[idx[1]] == 1
        #                 matrix = matrix + pol[k,1]*SetVars["MomMat"][1][1,1+idx[1]]*E;
        #             elseif length(idx) == 1 && v[idx[1]] == 2
        #                 matrix = matrix + pol[k,1]*SetVars["MomMat"][1][1+idx[1],1+idx[1]]*E;
        #             elseif length(idx) == 2 && v[idx[1]]+v[idx[2]] == 2
        #                 matrix = matrix + pol[k,1]*SetVars["MomMat"][1][1+idx[1],1+idx[2]]*E;
        #             else
        #                 for t = 1:length(SetVars["MomVar"])
        #                     for z = 1:length(SetVars["MomVar"][t])
        #                         if isequal(SetVars["MomSupp"][t][z, :], v)
        #                             # @printf("%d %d %d\n", i,j,z)
        #                             # add_to_expression!(matrix, pol[k,1]*vars["var"][z]*E)
        #                             matrix = matrix + pol[k,1]*SetVars["MomVar"][t][z]*E;
        #                             new = 0;
        #                             break
        #                         end
        #                     end
        #                     if new == 0
        #                         break
        #                     end
        #                 end
        #                 if new == 1
        #                     SetVars["var"] = vcat(SetVars["var"], @variable(model, [1:1, 1:1]));
        #                     SetVars["supp"] = vcat(SetVars["supp"], sparse(v)');
        #                     matrix = matrix + pol[k,1]*SetVars["var"][length(SetVars["var"])]*E;
        #                 end
        #             end
        #         end
        #     end
        # end
        # println(matrix)
    return matrix
end
# start = time();
# n = 2; order = 1; num = binomial(n+2*order+2, 2*order+2);
# basis = [1 0];
# pol = [1 1 0; -1 2 0];
# SetVars = Dict();
# model = Model(with_optimizer(Mosek.Optimizer));
# SetVars["var"] = @variable(model, y[1:num, 1:1]);
# SetVars["supp"] = gen_basis(n, 2*order+2);
# vars, matrix = LocalizationMatrix(model, pol, basis, order, SetVars);
# elapsed = time() - start;
