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
    OptVal, running_time, stat = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    return OptVal, running_time, stat
end
# @polyvar x[1:2];
# obj = x[1]*x[2];
# var = x; typ = "max";
# MomConst = Array{Any}(undef, 1);
# MomConst[1] = Dict(); MomConst[1]["basis"] = x; MomConst[1]["ord"] = 2;
# LocConst = Array{Any}(undef, 1);
# LocConst[1] = Dict(); LocConst[1]["basis"] = x; LocConst[1]["pol"] = 1-sum(x.^2); LocConst[1]["typ"] = ">="; LocConst[1]["ord"] = 1;
# options = Dict(); options["silent"] = true; options["quad"] = false;
# OptVal, running_time, status = solve_moment_auto(typ, var, obj, MomConst, LocConst, options);

function solve_moment_manual(typ, obj, MomConst, LocConst, options)
    running_time = Dict();
    start = time();
    model = Model(with_optimizer(Mosek.Optimizer));
    # model = Model(with_optimizer(DualOptimizer, Mosek.Optimizer()));
    NumMom = length(MomConst); NumLoc = length(LocConst); DimVar = size(obj, 2)-1;
    SetVars = Dict();
    SetVars["var"] = []#Array{GenericAffExpr{Float64,VariableRef},2}(undef, 1, 1);
    # SetVars["var"][:] = @variable(model, x[1:size(obj, 1), 1:1]);
    SetVars["supp"] = []#Array{Int64,2}(undef, 1, DimVar)#obj[:, 2:DimVar+1];
    SetVars["MomMat"] = Array{Any}(undef, NumMom); SetVars["MomConst"] = Array{Any}(undef, NumMom); SetVars["MomVar"] = Array{Any}(undef, NumMom); SetVars["MomSupp"] = Array{Any}(undef, NumMom);
    for i = 1:NumMom
        SetVars["MomMat"][i] = []; SetVars["MomConst"][i] = [];
    end
    SetVars["LocMat"] = Array{Any}(undef, NumLoc); SetVars["LocConst"] = Array{Any}(undef, NumLoc); SetVars["LocPol"] = Array{Any}(undef, NumLoc); #SetVars["LocVar"] = Array{Any}(undef, NumLoc);
    for i = 1:NumLoc
        de = size(LocConst[i]["pol"], 1)
        SetVars["LocMat"][i] = Array{Any}(undef, de);
        for j = 1:de
            SetVars["LocMat"][i][j] = []
        end
        SetVars["LocConst"][i] = []; SetVars["LocPol"][i] = [];
    end
    # turn = "off";
    # for i = 1:size(obj, 1)
    #     if isequal(obj[i, 2:DimVar+1]', zeros(1, DimVar))
    #         SetVars["var"][i] = 1;
    #         # @constraints(model, begin
    #         #     SetVars["var"][i] .== 1
    #         # end)
    #         turn = "on";
    #         break
    #     end
    # end
    if isequal(options["quad"], false)
    # if maximum(sum(obj[:, 2:DimVar+1], dims = 2)) > 2
        ObjQ = "no";
        objective = obj[:, 1]'*SetVars["var"][1:size(obj, 1), 1];
        if isequal(typ, "max")
            @objective(model, Max, objective);
        elseif isequal(typ, "min")
            @objective(model, Min, objective);
        end
    end
    for i = 1:NumMom
        if i == 1 && MomConst[i]["ord"] != 0
            # println(MomConst)
            MomMat = MomentMatrix(model, MomConst[i]["basis"], MomConst[i]["ord"], SetVars, ObjQuad = options["quad"]);
            SetVars["MomMat"][i] = MomMat; SetVars["MomConst"][i] = MomConst[i]["basis"]; SetVars["MomVar"][i] = SetVars["var"]; SetVars["MomSupp"][i] = SetVars["supp"];
            # println(MomMat)
            # SetVars["var"] = Vars["var"]; SetVars["supp"] = Vars["supp"];
            @printf("Building up moment matrices: %.2f%% (%d/%d)", i/NumMom*100, i, NumMom);
            if isequal(options["quad"], true)
                @constraint(model, MomMat[1,1] == 1)
                # println(size(obj)); println(size(SetVars["Mat"]))
                objective = sum(obj.*SetVars["MomMat"][i]);
                if isequal(typ, "max")
                    @objective(model, Max, objective);
                elseif isequal(typ, "min")
                    @objective(model, Min, objective);
                end
            end
            @constraint(model, Matrix(MomMat) in PSDCone())
        elseif MomConst[i]["ord"] != 0
            # println(MomConst)
            MomMat = MomentMatrix(model, MomConst[i]["basis"], MomConst[i]["ord"], SetVars, ObjQuad = false);
            SetVars["MomMat"][i] = MomMat; SetVars["MomConst"][i] = MomConst[i]["basis"]; SetVars["MomVar"][i] = SetVars["var"]; SetVars["MomSupp"][i] = SetVars["supp"];
            # println(MomMat)
            # SetVars["var"] = Vars["var"]; SetVars["supp"] = Vars["supp"];
            s1 = string(floor(Int64, (i-1)/NumMom*100)); s2 = string(i-1); s3 = string(NumMom);
            backstr = repeat("\b", 8+length(string(s1, s2, s3)));
            @printf("%s%.2f%% (%d/%d)", backstr, i/NumMom*100, i, NumMom);
            @constraint(model, Matrix(MomMat) in PSDCone())
        end
        # @printf("Building up moment matrices: %.2f%% (%d/%d)\n", i/NumMom*100, i, NumMom)
        # if isequal(turn, "off")
        #     for j = 1:size(SetVars["supp"], 1)
        #         if isequal(SetVars["supp"][j,:]', zeros(1, DimVar))
        #             SetVars["var"][j] = 1;
        #             # @constraints(model, begin
        #             #     SetVars["var"][j] == 1
        #             # end)
        #             turn = "on";
        #             break
        #         end
        #     end
        # end
        # MomMat = [];
    end
    @printf("\n")
    for i = 1:NumLoc
        # i = i+1
        LocMat = LocalizationMatrix(model, LocConst[i]["pol"], LocConst[i]["basis"], LocConst[i]["ord"], SetVars);
        SetVars["LocConst"][i] = LocConst[i]["basis"]; SetVars["LocPol"][i] = LocConst[i]["pol"];
        # println(LocMat)
        for j = 1:length(SetVars["TempLocMat"])
            SetVars["LocMat"][i][j] = SetVars["TempLocMat"][j];
        end
        # if i == 12
        #     println("")
        #     # println(SetVars["LocMat"][i][2])
        #     # println(vars["LocMat"][i][3][2,2]==vars["MomMat"][12][6,4])
        #     println(SetVars["MomMat"][3] .== SetVars["MomMat"][12])
        # end
        # println(LocConst[i]["pol"][:,1]); println(LocMat[1,1]);print("\n")
        # SetVars["var"] = Vars["var"]; SetVars["supp"] = Vars["supp"];
        if i == 1
            @printf("Building up localization matrices: %.2f%% (%d/%d)", i/NumMom*100, i, NumMom);
        else
            s1 = string(floor(Int64, (i-1)/NumLoc*100)); s2 = string(i-1); s3 = string(NumLoc);
            backstr = repeat("\b", 8+length(string(s1, s2, s3)));
            @printf("%s%.2f%% (%d/%d)", backstr, i/NumLoc*100, i, NumLoc);
        end
        # @printf("Building up localization matrices: %.2f%% (%d/%d)\n", i/NumLoc*100, i, NumLoc)
        if isequal(LocConst[i]["typ"], ">=")
            if size(LocMat, 1) == 1
                @constraint(model, LocMat .>= 0)
            elseif size(LocMat, 1) > 1
                @constraint(model, Matrix(LocMat) in PSDCone())
                # @constraint(model, LocMat[1,1] >= 0)
            end
        elseif isequal(LocConst[i]["typ"], "<=")
            if size(LocMat, 1) == 1
                @constraint(model, LocMat[1,1] .<= 0)
            elseif size(LocMat, 1) > 1
                @constraint(model, -Matrix(LocMat) in PSDCone())
                # @constraint(model, LocMat[1,1] <= 0)
            end
        elseif isequal(LocConst[i]["typ"], "==")
            @constraint(model, LocMat .== 0)
            # @constraint(model, LocMat[1,1] == 0)
        end
        # if isequal(turn, "off")
        #     for j = 1:size(SetVars["supp"], 1)
        #         if isequal(SetVars["supp"][j,:]', zeros(1, DimVar))
        #             @constraints(model, begin
        #                 SetVars["var"][j] == 1
        #             end)
        #             println("Aha")
        #             turn = "on";
        #             break
        #         end
        #     end
        # end
        # LocMat = [];
    end
    @printf("\n")
    # println(objective)
    running_time["model"] = time() - start;
    @printf("Problem type: %s\n", uppercase(typ))
    @printf("%d moment matrices, %d localization matrices, %d scalar variables, %d constraints\n", NumMom, NumLoc, length(SetVars["var"]), NumMom+NumLoc); #SetVars = [];
    @printf("Solver started.............solving.............")
    MOI.set(model, MOI.Silent(), options["silent"]);
    optimize!(model); #println(value.(SetVars["MomMat"][1][1,:]))
    @printf("Solver finished.\n")
    running_time["solv"] = solve_time(model);
    OptVal = objective_value(model);
    status = termination_status(model);
    @printf("Solution is: %.2f, ", OptVal)
    @printf("solver status is: "); println(status)
    # if status == MOI.OPTIMAL
    #     @printf("The problem is successfully solved! Optimal value is: %.2f\n", OptVal)
    # else
    #     @printf("The solver encountered some issues: ")
    #     println(termination_status(model))
    # end
    @printf("Total running time is: %.2f seconds (modeling time: %.2f seconds, solving time: %.2f seconds)\n", running_time["model"] + running_time["solv"], running_time["model"], running_time["solv"])
    # println(value.(SetVars["var"]))
    return OptVal, running_time, status, value.(SetVars["MomMat"][1])
end
# typ = "max"; obj = [1 1 1];
# options = Dict();
# options["silent"] = true; options["quad"] = false;
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
