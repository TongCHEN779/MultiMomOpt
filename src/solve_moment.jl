function solve_moment_manual(typ, obj, MomConst, LocConst, options; cliques = [], sizes = [])
    running_time = Dict();
    start = time();
    # for older JuMP version: model = Model(with_optimizer(Mosek.Optimizer))
    model = Model(Mosek.Optimizer);
    NumMom = length(MomConst); # number of moment matrices
    NumLoc = length(LocConst); # number of localizing matrices
    # store the moment and localizing matrices
    vars = Dict();
    vars["MomMat"] = Array{Any}(undef, NumMom);
    vars["LocMat"] = Array{Any}(undef, NumLoc);
    # generate moment matrices
    @printf("Building up moment matrices: \n");
    for i in ProgressBar(1:NumMom+1)
        if i == NumMom+1
            continue # to avoid the bug of ProgressBar
        end
        MomMat = MomentMatrix(model, MomConst[i]["basis"], MomConst[i]["ord"], MomConst[1:i-1], vars; ObjQuad = options["quad"]);
        vars["MomMat"][i] = MomMat;
        if i == 1
            @constraints(model, begin
                Matrix(MomMat) in PSDCone()
                MomMat[1,1] == 1
            end)
        else
            @constraint(model, Matrix(MomMat) in PSDCone())
        end
    end
    # generate localizing matrices
    @printf("Building up localizing matrices: \n");
    for i in ProgressBar(1:NumLoc)
        LocMat, LocMatSep = LocalizingMatrix(model, LocConst[i]["pol"], LocConst[i]["basis"], LocConst[i]["ord"], MomConst, LocConst[1:i-1], vars);
        vars["LocMat"][i] = LocMatSep;
        if isequal(LocConst[i]["typ"], ">=") && size(LocMat, 1) == 1
            @constraint(model, LocMat .>= 0)
        elseif isequal(LocConst[i]["typ"], ">=") && size(LocMat, 1) > 1
            @constraint(model, Matrix(LocMat) in PSDCone())
        elseif isequal(LocConst[i]["typ"], "<=") && size(LocMat, 1) == 1
            @constraint(model, LocMat[1,1] .<= 0)
        elseif isequal(LocConst[i]["typ"], "<=") && size(LocMat, 1) > 1
            @constraint(model, -Matrix(LocMat) in PSDCone())
        elseif isequal(LocConst[i]["typ"], "==")
            @constraint(model, LocMat .== 0)
        end
    end
    # generate objective polynomial
    if isequal(options["quad"], true) && isequal(options["clique"], "on") # sparse quadratic
        n = size(obj,1)
        objective = 0
        for i = 1:n
            for j = i:n
                if obj[i,j] != 0
                    for k = 1:length(sizes)
                        Clq = vcat(1, cliques[k]'.+1)
                        if i == j && i in Clq
                            idx = findall(Clq .== i)
                            objective += obj[i,i]*vars["MomMat"][k][idx[1], idx[1]]
                            break
                        elseif i != j && i in Clq && j in Clq
                            idx1 = findall(Clq .== i)
                            idx2 = findall(Clq .== j)
                            objective += 2*obj[i,j]*vars["MomMat"][k][idx1[1], idx2[1]]
                            break
                        end
                    end
                end
            end
        end
        if isequal(typ, "max")
            @objective(model, Max, objective);
        elseif isequal(typ, "min")
            @objective(model, Min, objective);
        end
    elseif isequal(options["quad"], true) && isequal(options["clique"], "off") # dense quadratic
        objective = sum(obj.*vars["MomMat"][1]);
        if isequal(typ, "max")
            @objective(model, Max, objective);
        elseif isequal(typ, "min")
            @objective(model, Min, objective);
        end
    else isequal(options["quad"], false) # non-quadratic optimization
        error("NotImplementedError")
    end
    # solve the problem
    @printf("\n")
    running_time["model"] = time() - start;
    @printf("Problem type: %s\n", uppercase(typ))
    @printf("%d moment matrices, %d localizing matrices.\n", NumMom, NumLoc);
    @printf("Solver started.............solving.............")
    MOI.set(model, MOI.Silent(), options["silent"]);
    optimize!(model);
    @printf("Solver finished.\n")
    running_time["solv"] = solve_time(model);
    OptVal = objective_value(model);
    status = termination_status(model);
    @printf("Solution is: %.2f, ", OptVal)
    @printf("solver status is: "); println(status)
    @printf("Total running time is: %.2f seconds (modeling time: %.2f seconds, solving time: %.2f seconds)\n", running_time["model"] + running_time["solv"], running_time["model"], running_time["solv"])
    if isequal(options["clique"], "on") # sparse
        M = Array{Any}(undef, length(sizes));
        for i = 1:length(sizes)
            M[i] = value.(vars["MomMat"][i])
        end
        return OptVal, running_time, status, M
    elseif isequal(options["clique"], "off") # dense
        return OptVal, running_time, status, value.(vars["MomMat"][1])
    end
end
