function solve_moment_maxcut(A, W, options)
    n = size(A, 1); typ = "max";
    L = sparse(diagm(((A.*W)*ones(n,1))[:,1])) - A.*W;
    E0 = sparse([zeros(1+n, n-1) [2 zeros(1,n); ones(n,1) I(n)]]);
    E = sparse([I(n); I(n)]);
    L = L/2;
    for i = 1:n
        L[i,i] = L[i,i]/2;
    end
    obj = sparse([zeros(n+1,1) E0[1:1+n, n+1:2*n]]);
    for i = 2:n+1
        obj = vcat(obj, [L[i-1:n,i-1] E0[1:n-i+2, n-i+2:2*n-i+1]]);
    end
    if options["level"] == 0 # Shor's relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        for i = 1:n
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[i] = Dict();
            LocConst[i]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
            LocConst[i]["basis"] = sparse(I(n));
            LocConst[i]["typ"] = "==";
            LocConst[i]["ord"] = 0;
        end
    elseif options["level"] > 0 && isequal(options["clique"], "on")
        cliques, sizes = chordal_extension(A);
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
        NumMom = 1;
        for i = 1:length(cliques)
            clique = [cliques[i] cliques[i]];
            if options["level"] < sizes[i]
                for j = 1:sizes[i]
                    NumMom = NumMom + 1;
                    MomConst = vcat(MomConst, Array{Any}(undef, 1));
                    MomConst[NumMom] = Dict();
                    MomConst[NumMom]["basis"] = sparse(E[clique[j:j+options["level"]-1], :]');
                    MomConst[NumMom]["ord"] = options["order"];
                end
            else
                NumMom = NumMom + 1;
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[cliques[i], :]');
                MomConst[NumMom]["ord"] = options["order"];
            end
        end
        LocConst = Array{Any}(undef, 0);
        NumLoc = 0;
        for i = 1:n
            for j = 1:length(cliques)
                if i in cliques[j]
                    clique = [cliques[j] cliques[j]]; idx = j;
                    break
                end
            end
            if options["level"] < sizes[idx]
                for j = 1:sizes[idx]
                    NumLoc = NumLoc + 1;
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                    LocConst[NumLoc]["basis"] = sparse(E[clique[j:j+options["level"]-1], :]');
                    LocConst[NumLoc]["typ"] = "==";
                    LocConst[NumLoc]["ord"] = options["order"] - 1;
                end
            else
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[cliques[idx], :]');
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = options["order"] - 1;
            end
        end
    elseif options["level"] > 0 && isequal(options["clique"], "off")
        MomConst = Array{Any}(undef, 1); LocConst = Array{Any}(undef, 0);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
        NumMom = 1; NumLoc = 0;
        clique = [1:n; 1:n]';
        if options["level"] < n
            for i = 1:n
                NumMom = NumMom + 1;
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[clique[i:i+options["level"]-1], :]');
                MomConst[NumMom]["ord"] = options["order"];
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[clique[i:i+options["level"]-1], :]');
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = options["order"] - 1;
            end
        else
            NumMom = NumMom + 1;
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[1:n, :]');
            MomConst[NumMom]["ord"] = options["order"];
            for i = 1:n
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[1:n, :]');
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = options["order"] - 1;
            end
        end
    end
    @printf("\nMAX-CUT problem: %d vertices, %d edges, duplicated %s, clique %s, order %d, level %d\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), uppercase(options["duplicated"]), uppercase(options["clique"]), options["order"], options["level"])
    # OptVal, time, stat = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    # return OptVal, time, stat
end
# vars = matread("A.mat");
# A = vars["A"]; W = ones(size(A, 1), size(A, 1));
# options = Dict();
# options["duplicated"] = "on"; options["level"] = 0; options["clique"] = "on"; options["ord"] = 2;
# solve_moment_maxcut(A, W, options)

function solve_moment_mip(A, options)
    n = size(A, 1); typ = "min";
    E0 = sparse([zeros(1+n, n-1) [2 zeros(1,n); ones(n,1) I(n)]]);
    E = sparse([I(n); I(n)]);
    A = 2*A;
    for i = 1:n
        A[i,i] = A[i,i]/2;
    end
    obj = sparse([zeros(n+1,1) E0[1:1+n, n+1:2*n]]);
    for i = 2:n+1
        obj = vcat(obj, [A[i-1:n,i-1] E0[1:n-i+2, n-i+2:2*n-i+1]]);
    end
    if options["level"] == 0 # Shor's relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        for i = 1:n
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[i] = Dict();
            LocConst[i]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
            LocConst[i]["basis"] = sparse(I(n));
            LocConst[i]["typ"] = "==";
            LocConst[i]["ord"] = 0;
        end
    elseif options["level"] > 0 && isequal(options["clique"], "on")
        cliques, sizes = chordal_extension(A);
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
        NumMom = 1;
        for i = 1:length(cliques)
            clique = [cliques[i] cliques[i]];
            if options["level"] < sizes[i]
                for j = 1:sizes[i]
                    NumMom = NumMom + 1;
                    MomConst = vcat(MomConst, Array{Any}(undef, 1));
                    MomConst[NumMom] = Dict();
                    MomConst[NumMom]["basis"] = sparse(E[clique[j:j+options["level"]-1], :]');
                    MomConst[NumMom]["ord"] = options["order"];
                end
            else
                NumMom = NumMom + 1;
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[cliques[i], :]');
                MomConst[NumMom]["ord"] = options["order"];
            end
        end
        LocConst = Array{Any}(undef, 0);
        NumLoc = 0;
        for i = 1:n
            for j = 1:length(cliques)
                if i in cliques[j]
                    clique = [cliques[j] cliques[j]]; idx = j;
                    break
                end
            end
            if options["level"] < sizes[idx]
                for j = 1:sizes[idx]
                    NumLoc = NumLoc + 1;
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                    LocConst[NumLoc]["basis"] = sparse(E[clique[j:j+options["level"]-1], :]');
                    LocConst[NumLoc]["typ"] = "==";
                    LocConst[NumLoc]["ord"] = options["order"] - 1;
                end
            else
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[cliques[idx], :]');
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = options["order"] - 1;
            end
        end
    elseif options["level"] > 0 && isequal(options["clique"], "off")
        MomConst = Array{Any}(undef, 1); LocConst = Array{Any}(undef, 0);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
        NumMom = 1; NumLoc = 0;
        clique = [1:n; 1:n]';
        if options["level"] < n
            for i = 1:n
                NumMom = NumMom + 1;
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[clique[i:i+options["level"]-1], :]');
                MomConst[NumMom]["ord"] = options["order"];
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[clique[i:i+options["level"]-1], :]');
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = options["order"] - 1;
            end
        else
            NumMom = NumMom + 1;
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[1:n, :]');
            MomConst[NumMom]["ord"] = options["order"];
            for i = 1:n
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[1:n, :]');
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = options["order"] - 1;
            end
        end
    end
    @printf("\nMixed Integer Programming (MIP): %d variables, %d entries, duplicated %s, clique %s, order %d, level %d\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), uppercase(options["duplicated"]), uppercase(options["clique"]), options["order"], options["level"])
    # OptVal, time, stat = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    # return OptVal, time, stat
end
# vars = matread("A.mat");
# A = vars["A"];
# options = Dict();
# options["duplicated"] = "on"; options["level"] = 0; options["clique"] = "on"; options["ord"] = 2;
# solve_moment_mip(A, options)

function solve_moment_auto(typ, var, obj, MomConst, LocConst, options)
    obj = pol2supp(obj, var);
    NumMom = length(MomConst); NumLoc = length(LocConst);
    for i = 1:NumMom
        MomConst[i]["basis"] = basis2supp(MomConst[i]["basis"], var);
    end
    for i = 1:NumLoc
        LocConst[i]["pol"] = pol2supp(LocConst[i]["pol"], var);
        LocConst[i]["basis"] = basis2supp(LocConst[i]["basis"], var);
    end
    OptVal, time, stat = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    return OptVal, time, stat
end
# @polyvar x[1:3];
# obj = 1.5 + x[1] + x[2]^2*x[3];
# var = x; typ = "min";
# MomConst = Array{Any}(undef, 1);
# MomConst[1] = Dict(); MomConst[1]["basis"] = x;
# LocConst = Array{Any}(undef, 1);
# LocConst[1] = Dict(); LocConst[1]["basis"] = x; LocConst[1]["pol"] = sum(x);
# options = Dict();

function solve_moment_manual(typ, obj, MomConst, LocConst, options)
    running_time = Dict();
    start = time();
    model = Model(with_optimizer(Mosek.Optimizer));
    NumMom = length(MomConst); NumLoc = length(LocConst); DimVar = size(obj, 2)-1;
    SetVars = Dict();
    SetVars["var"] = @variable(model, [1:size(obj, 1), 1:1]);
    SetVars["supp"] = obj[:, 2:DimVar+1];
    turn = "off";
    for i = 1:size(obj, 1)
        if isequal(obj[i, 2:DimVar+1]', zeros(1, DimVar))
            @constraints(model, begin
                SetVars["var"][i] == 1
            end)
            turn = "on";
            break
        end
    end
    if maximum(sum(obj[:, 2:DimVar+1], dims = 2)) > 2
        ObjQ = "no";
        objective = obj[:, 1]'*SetVars["var"][1:size(obj, 1), 1];
        if isequal(typ, "max")
            @objective(model, Max, objective);
        elseif isequal(typ, "min")
            @objective(model, Min, objective);
        end
    else
        ObjQ = "yes";
    end
    for i = 1:NumMom
        @printf("Building up moment matrices: %.2f%% (%d/%d)\n", i/NumMom*100, i, NumMom)
        if i == 1
            # @printf("Building up moment matrices: %.2f%% (%d/%d)\n", i/NumMom*100, i, NumMom)
            SetVars, MomMat = MomentMatrix(model, MomConst[i]["basis"], MomConst[i]["ord"], SetVars, ObjQuad = ObjQ);
            if isequal(ObjQ, "yes")
                objective = obj[:, 1]'*SetVars["var"][1:size(obj, 1), 1];
                if isequal(typ, "max")
                    @objective(model, Max, objective);
                elseif isequal(typ, "min")
                    @objective(model, Min, objective);
                end
            end
        else
            # s1 = string(floor(Int64, (i-1)/NumMom*100)); s2 = string(i-1); s3 = string(NumMom);
            # backstr = string(repeat("\b", 9+length([s1,s2,s3])));
            # @printf("%s%.2f%% (%d/%d)\n", backstr, i/NumMom*100, i, NumMom)
            SetVars, MomMat = MomentMatrix(model, MomConst[i]["basis"], MomConst[i]["ord"], SetVars, ObjQuad = "no");
        end
        @constraints(model, begin
            MomMat in PSDCone()
        end)
        if isequal(turn, "off")
            for j = 1:size(SetVars["supp"], 1)
                if isequal(SetVars["supp"][j,:]', zeros(1, DimVar))
                    @constraints(model, begin
                        SetVars["var"][j] == 1
                    end)
                    turn = "on";
                    break
                end
            end
        end
    end
    for i = 1:NumLoc
        @printf("Building up localization matrices: %.2f%% (%d/%d)\n", i/NumLoc*100, i, NumLoc)
        # if i == 1
        #     fprintf('Building up localization matrices: %.2f%% (%d/%d)\n', i/NumLoc*100, i, NumLoc)
        # else
        #     s1 = num2str(floor((i-1)/NumLoc*100)); s2 = num2str(i-1); s3 = num2str(NumLoc);
        #     fprintf([repmat('\b', 1, 9+length([s1,s2,s3])), '%.2f%% (%d/%d)\n'], i/NumLoc*100, i, NumLoc)
        # end
        SetVars, LocMat = LocalizationMatrix(model, LocConst[i]["pol"], LocConst[i]["basis"], LocConst[i]["ord"], SetVars);
        if isequal(LocConst[i]["typ"], ">=")
            @constraints(model, begin
                LocMat in PSDCone()
            end)
        elseif isequal(LocConst[i]["typ"], "<=")
            @constraints(model, begin
                -LocMat in PSDCone()
            end)
        elseif isequal(LocConst[i]["typ"], "==")
            @constraints(model, begin
                -LocMat .== 0
            end)
        end
        if isequal(turn, "off")
            for j = 1:size(SetVars["supp"], 1)
                if isequal(SetVars["supp"][j,:]', zeros(1, DimVar))
                    @constraints(model, begin
                        SetVars["var"][j] == 1
                    end)
                    turn = "on";
                    break
                end
            end
        end
    end
    running_time["model"] = time() - start;
    @printf("Problem type: %s\n", uppercase(typ))
    @printf("%d moment matrices, %d localization matrices, %d variables, %d constraints\n", NumMom, NumLoc, length(SetVars["var"]), NumMom+NumLoc)
    @printf("Solver started.............solving.............")
    MOI.set(model, MOI.Silent(), options["silent"]);
    optimize!(model);
    @printf("Solver finished.\n")
    running_time["solv"] = solve_time(model);
    OptVal = objective_value(model);
    status = termination_status(model);
    if status == MOI.OPTIMAL
        @printf("The problem is successfully solved! Optimal value is: %.2f\n", OptVal)
    else
        @printf("The solver encountered some issues: ")
        println(termination_status(model))
    end
    @printf("Total running time is: %.2f seconds (modeling time: %.2f seconds, solver time: %.2f seconds)\n", running_time["model"] + running_time["solv"], running_time["model"], running_time["solv"])
    return OptVal, running_time, status
end
# typ = "min"; obj = [0 0 0; 0 1 0; 0 0 1; 0 2 0; 1 1 1; 0 0 2];
# options = Dict();
# options["level"] = 0; options["clique"] = "on"; options["silent"] = true;
# MomConst = Array{Any}(undef, 1);
# MomConst[1] = Dict();
# MomConst[1]["basis"] = [1 0; 0 1];
# MomConst[1]["ord"] = 2;
# LocConst = Array{Any}(undef, 1);
# LocConst[1] = Dict();
# LocConst[1]["basis"] = [1 0; 0 1]; 
# LocConst[1]["pol"] = [1 0 0; -1 2 0; -1 0 2];
# LocConst[1]["typ"] = ">=";
# LocConst[1]["ord"] = 1;
# OptVal, running_time, status = solve_moment_manual(typ, obj, MomConst, LocConst, options);
