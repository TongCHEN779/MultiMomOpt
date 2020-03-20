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
    NumMom = length(MomConst); NumLoc = length(LocConst); DimVar = size(obj, 2)-1;
    SetVars = Dict();
    SetVars["var"] = Array{GenericAffExpr{Float64,VariableRef},2}(undef, size(obj, 1), 1);
    SetVars["var"][:] = @variable(model, x[1:size(obj, 1), 1:1]);
    SetVars["supp"] = obj[:, 2:DimVar+1];
    turn = "off";
    for i = 1:size(obj, 1)
        if isequal(obj[i, 2:DimVar+1]', zeros(1, DimVar))
            SetVars["var"][i] = 1;
            # @constraints(model, begin
            #     SetVars["var"][i] .== 1
            # end)
            turn = "on";
            break
        end
    end
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
        if i == 1
            SetVars, MomMat = MomentMatrix(model, MomConst[i]["basis"], MomConst[i]["ord"], SetVars, ObjQuad = options["quad"]); MomMat[1,1] = 1;
            @printf("Building up moment matrices: %.2f%% (%d/%d)", i/NumMom*100, i, NumMom);
            if isequal(options["quad"], true)
                objective = obj[:, 1]'*SetVars["var"][1:size(obj, 1), 1];
                if isequal(typ, "max")
                    @objective(model, Max, objective);
                elseif isequal(typ, "min")
                    @objective(model, Min, objective);
                end
            end
        else
            SetVars, MomMat = MomentMatrix(model, MomConst[i]["basis"], MomConst[i]["ord"], SetVars, ObjQuad = false); MomMat[1,1] = 1;
            s1 = string(floor(Int64, (i-1)/NumMom*100)); s2 = string(i-1); s3 = string(NumMom);
            backstr = repeat("\b", 8+length(string(s1, s2, s3)));
            @printf("%s%.2f%% (%d/%d)", backstr, i/NumMom*100, i, NumMom);
        end
        # @printf("Building up moment matrices: %.2f%% (%d/%d)\n", i/NumMom*100, i, NumMom)
        @constraint(model, MomMat in PSDCone())
        if isequal(turn, "off")
            for j = 1:size(SetVars["supp"], 1)
                if isequal(SetVars["supp"][j,:]', zeros(1, DimVar))
                    SetVars["var"][j] = 1;
                    # @constraints(model, begin
                    #     SetVars["var"][j] == 1
                    # end)
                    turn = "on";
                    break
                end
            end
        end
    end
    @printf("\n")
    for i = 1:NumLoc
        SetVars, LocMat = LocalizationMatrix(model, LocConst[i]["pol"], LocConst[i]["basis"], LocConst[i]["ord"], SetVars);
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
                @constraint(model, LocMat in PSDCone())
            end
        elseif isequal(LocConst[i]["typ"], "<=")
            if size(LocMat, 1) == 1
                @constraint(model, LocMat .<= 0)
            elseif size(LocMat, 1) > 1
                @constraint(model, -LocMat in PSDCone())
            end
        elseif isequal(LocConst[i]["typ"], "==")
            @constraint(model, LocMat .== 0)
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
    end
    @printf("\n")
    # println(objective)
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
    @printf("Solution is: %.2f, ", OptVal)
    @printf("solver status is: "); println(termination_status(model))
    # if status == MOI.OPTIMAL
    #     @printf("The problem is successfully solved! Optimal value is: %.2f\n", OptVal)
    # else
    #     @printf("The solver encountered some issues: ")
    #     println(termination_status(model))
    # end
    @printf("Total running time is: %.2f seconds (modeling time: %.2f seconds, solving time: %.2f seconds)\n", running_time["model"] + running_time["solv"], running_time["model"], running_time["solv"])
    # println(value.(SetVars["var"]))
    return OptVal, running_time, status
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
