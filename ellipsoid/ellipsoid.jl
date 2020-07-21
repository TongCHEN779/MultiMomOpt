using Plots
using SumOfSquares
using DynamicPolynomials
using MosekTools
using LinearAlgebra
using Printf
using MathOptInterface

# compute the outer approximation ellipsoid
function OuterApproximation(Q0, b0, c0, ord, P, ξ, method)
    p1 = size(P, 2); p2 = size(P, 1)
    model = SOSModel(with_optimizer(Mosek.Optimizer))
    # objective variable
    Q = @variable(model, [1:p2, 1:p2], Symmetric);
    b = @variable(model, [1:p2, 1:1]);
    c = @variable(model, [1:1, 1:1]);
    if isequal(method, "tong")
        @polyvar x[1:p1] z[1:p2] lower_bound
        @variable(model, lower_bound)
        @objective(model, Max, lower_bound)
        Mat = [-Q b/2; b'/2 1-c[1]];
        vect = lower_bound
        for i = 1:p2
            vect = hcat(vect, -Q[i,i:p2]')
        end
        @constraints(model, begin
            Mat in PSDCone()
            vect[:] in MOI.RootDetConeTriangle(p2)
        end)
        X0 = monomials([x; z], 0:2*ord); σ0 = @variable(model, [1:1], Poly(X0))
        X1 = monomials([x; z], 0:2*(ord-1)); σ1 = @variable(model, [1:1], Poly(X1))
        obj = z'*Q*z.+b'*z.+c - σ0[1] - σ1[1] * (x'*Q0*x.+b0'*x.+c0)
        Λ = @variable(model, [1:p2,1:p2], Symmetric)
        for i = 1:p2-1
            for j = i+1:p2
                obj = obj + Λ[i,j]*((z[j] - z[i])^2 - (z[j] - z[i])*(P[j,:]'*x + ξ[j] - P[i,:]'*x - ξ[i]))
            end
        end
        @constraints(model, begin
            σ0 .>= 0; σ1 .>= 0; Λ .>= 0
        end)
        for i = 1:p2
            X2 = monomials([x; z], 0:2*(ord-1)); σ2 = @variable(model, [1:1], Poly(X2))
            X3 = monomials([x; z], 0:2*(ord-1)); σ3 = @variable(model, [1:1], Poly(X3))
            X4 = monomials([x; z], 0:2*(ord-1)); σ4 = @variable(model, [1:1], Poly(X4))
            obj = obj - σ2[1] * (z[i]*(z[i]-P[i,:]'*x-ξ[i])) - σ3[1] * z[i] - σ4[1] * (z[i]-P[i,:]'*x-ξ[i])
            @constraints(model, begin
                σ3 .>= 0; σ4 .>= 0;
            end)
        end
        # println("Shor:")
        # println((-obj-σ0[1])[1])
        @constraint(model, obj .== 0)
    elseif isequal(method, "split")
        @polyvar x[1:p2] z[1:p2] lower_bound
        @variable(model, lower_bound)
        @objective(model, Max, lower_bound)
        Mat = [-Q b/2; b'/2 1-c[1]];
        vect = lower_bound
        for i = 1:p2
            vect = hcat(vect, -Q[i,i:p2]')
        end
        @constraints(model, begin
            Mat in PSDCone()
            vect[:] in MOI.RootDetConeTriangle(p2)
        end)
        X0 = monomials([x; z], 0:2*ord); σ0 = @variable(model, [1:1], Poly(X0))
        X1 = monomials([x; z], 0:2*(ord-1)); σ1 = @variable(model, [1:1], Poly(X1))
        D = svd(P); Q0 = D.U*diagm(D.S)^(-1)*D.V'*Q0*D.V*diagm(D.S)^(-1)*D.U'; b0 = D.U*diagm(D.S)^(-1)*D.V'*b0; c0 = c0
        obj = z'*Q*z.+b'*z.+c - σ0[1] - σ1[1] * ((x-ξ)'*Q0*(x-ξ).+b0'*(x-ξ).+c0)
        Λ = @variable(model, [1:p2,1:p2], Symmetric)
        for i = 1:p2-1
            for j = i+1:p2
                obj = obj + Λ[i,j]*((z[j] - z[i])^2 - (z[j] - z[i])*(x[j] - x[i]))
            end
        end
        @constraints(model, begin
            σ0 .>= 0; σ1 .>= 0; Λ .>= 0
        end)
        for i = 1:p2
            X2 = monomials([x; z], 0:2*(ord-1)); σ2 = @variable(model, [1:1], Poly(X2))
            X3 = monomials([x; z], 0:2*(ord-1)); σ3 = @variable(model, [1:1], Poly(X3))
            X4 = monomials([x; z], 0:2*(ord-1)); σ4 = @variable(model, [1:1], Poly(X4))
            obj = obj - σ2[1] * (z[i]*(z[i]-x[i])) - σ3[1] * z[i] - σ4[1] * (z[i]-x[i])
            @constraints(model, begin
                σ3 .>= 0; σ4 .>= 0;
            end)
        end
        @constraint(model, obj .== 0)
    elseif isequal(method, "jean")
        @polyvar x[1:p1] z[1:p2] lower_bound
        @variable(model, lower_bound)
        @objective(model, Max, lower_bound)
        Mat = [-Q b/2; b'/2 1-c[1]];
        vect = lower_bound
        for i = 1:p2
            vect = hcat(vect, -Q[i,i:p2]')
        end
        @constraints(model, begin
            Mat in PSDCone()
            vect[:] in MOI.RootDetConeTriangle(p2)
        end)
        X0 = monomials([x; z], 0:2*ord); σ0 = @variable(model, [1:1], Poly(X0))
        X1 = monomials([x; z], 0:2*(ord-1)); σ1 = @variable(model, [1:1], Poly(X1))
        obj = 1/4*(z+P*x+ξ)'*Q*(z+P*x+ξ).+1/2*b'*(z+P*x+ξ).+c - σ0[1] - σ1[1] * (x'*Q0*x.+b0'*x.+c0)
        @constraints(model, begin
            σ0 .>= 0; σ1 .>= 0;
        end)
        for i = 1:p2
            X2 = monomials([x; z], 0:2*(ord-1)); σ2 = @variable(model, [1:1], Poly(X2))
            X3 = monomials([x; z], 0:2*(ord-1)); σ3 = @variable(model, [1:1], Poly(X3))
            obj = obj - σ2[1] * (z[i]^2-(P[i,:]'*x+ξ[i])^2) - σ3[1] * z[i]
            @constraints(model, begin
                σ3 .>= 0;
            end)
        end
        @constraint(model, obj .== 0)
    elseif isequal(method, "hr-2")
        @polyvar x[1:p2] z[1:p2] lower_bound
        @variable(model, lower_bound)
        @objective(model, Max, lower_bound)
        Mat = [-Q b/2; b'/2 1-c[1]];
        vect = lower_bound
        for i = 1:p2
            vect = hcat(vect, -Q[i,i:p2]')
        end
        @constraints(model, begin
            Mat in PSDCone()
            vect[:] in MOI.RootDetConeTriangle(p2)
        end)
        X0 = monomials([x; z], 0:2*1); σ0 = @variable(model, [1:1], Poly(X0))
        X01 = monomials(x, 0:2*ord); σ01 = @variable(model, [1:1], Poly(X01))
        X1 = monomials(x, 0:2*(ord-1)); σ1 = @variable(model, [1:1], Poly(X1))
        obj = z'*Q*z.+b'*z.+c - σ0[1] - σ01[1] - σ1[1] * (x'*Q0*x.+b0'*x.+c0)
        Λ = @variable(model, [1:p2,1:p2], Symmetric)
        for i = 1:p2-1
            for j = i+1:p2
                if i != j
                    obj = obj + Λ[i,j]*((z[j] - z[i])^2 - (z[j] - z[i])*(P[j,:]'*x + ξ[j] - P[i,:]'*x - ξ[i]))
                end
            end
        end
        @constraints(model, begin
            σ0 .>= 0; σ01 .>= 0; σ1 .>= 0; Λ .>= 0
        end)
        for i = 1:p2
            X02 = monomials([x[i]; z[i]], 0:2*ord); σ02 = @variable(model, [1:1], Poly(X02))
            X2 = monomials([x[i]; z[i]], 0:2*(ord-1)); σ2 = @variable(model, [1:1], Poly(X2))
            X3 = monomials([x[i]; z[i]], 0:2*(ord-1)); σ3 = @variable(model, [1:1], Poly(X3))
            X4 = monomials([x[i]; z[i]], 0:2*(ord-1)); σ4 = @variable(model, [1:1], Poly(X4))
            obj = obj - σ02[1] - σ2[1] * (z[i]*(z[i]-P[i,:]'*x-ξ[i])) - σ3[1] * z[i] - σ4[1] * (z[i]-P[i,:]'*x-ξ[i])
            @constraints(model, begin
                σ02 .>= 0; σ3 .>= 0; σ4 .>= 0;
            end)
        end
        @constraint(model, obj .== 0)
    end
    # solve
    MOI.set(model, MOI.Silent(), true);
    optimize!(model)
    @printf("Optimal value: %.2f, termination status: %s, primal status: %s, dual status: %s.\n", value.(lower_bound), termination_status(model), primal_status(model), dual_status(model))
    # println(value.(Q)); println(value.(b)); println(value.(c))
    # return value.(-Q*Q), value.(-2 .* Q*b), value.(1 .- b'*b)
    return value.(Q), value.(b), value.(c), value.(lower_bound)
end

function OuterApproximationSublevel(Q0, b0, c0, P, ξ, lv; method="cycle")
    l = length(P);
    p = Int.(zeros(l+1,1)); p[1] = size(P[1],2);
    for i = 1:l
        p[i+1] = size(P[i],1)
    end
    dx = Int.(zeros(l+2,1));
    for i = 1:l+1
        dx[i+1] = sum(p[1:i]);
    end
    model = SOSModel(with_optimizer(Mosek.Optimizer))
    # objective variable
    Q = @variable(model, [1:p[l+1], 1:p[l+1]], Symmetric);
    b = @variable(model, [1:p[l+1], 1:1]);
    c = @variable(model, [1:1, 1:1]);
    @polyvar lower_bound
    @variable(model, lower_bound)
    @objective(model, Max, lower_bound)
    Mat = [-Q b/2; b'/2 1-c[1]];
    vect = lower_bound
    for i = 1:p[l+1]
        vect = hcat(vect, -Q[i,i:p[l+1]]')
    end
    @constraints(model, begin
        Mat in PSDCone()
        vect[:] in MOI.RootDetConeTriangle(p[l+1])
    end)
    @polyvar x[1:dx[l+2]]
    X0 = monomials(x, 0:2*1); σ0 = @variable(model, [1:1], Poly(X0))
    obj = x[dx[l+1]+1:dx[l+2]]'*Q*x[dx[l+1]+1:dx[l+2]].+b'*x[dx[l+1]+1:dx[l+2]].+c - σ0[1]
    @constraint(model, σ0 .>= 0)
    if method == "cycle"
        for k = dx[1]:dx[2]-lv
            # X00 = monomials(x[k+1:k+lv], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00))
            X1 = monomials(x[k+1:k+lv], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1))
            obj = obj- σ1[1] * (x[1:dx[2]]'*Q0*x[1:dx[2]].+b0'*x[1:dx[2]].+c0)
            @constraint(model, σ1 .>= 0)
        end
        for k = dx[2]-lv+1:dx[2]-1
            # X00 = monomials(x[k+1:k+lv], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00))
            X1 = monomials(x[1:lv-dx[2]+k], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1))
            obj = obj- σ1[1] * (x[1:dx[2]]'*Q0*x[1:dx[2]].+b0'*x[1:dx[2]].+c0)
            @constraint(model, σ1 .>= 0)
        end
        for i = 1:l
            for j = 1:p[i+1]
                for k = dx[i]:dx[i+1]-lv
                    X00 = monomials([x[k+1:k+lv]; x[dx[i+1]+j]], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00));
                    X1 = monomials([x[k+1:k+lv]; x[dx[i+1]+j]], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1));
                    X2 = monomials([x[k+1:k+lv]; x[dx[i+1]+j]], 0:2*1); σ2 = @variable(model, [1:1], Poly(X2));
                    X3 = monomials([x[k+1:k+lv]; x[dx[i+1]+j]], 0:2*1); σ3 = @variable(model, [1:1], Poly(X3));
                    obj = obj - σ00[1] - σ1[1] * (x[dx[i+1]+j]*(x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])) - σ2[1] * x[dx[i+1]+j] - σ3[1] * (x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])
                    @constraints(model, begin
                        σ00 .>= 0; σ2 .>= 0; σ3 .>= 0;
                    end)
                end
                for k = dx[i+1]-lv+1:dx[i+1]-1
                    X00 = monomials([x[1:lv-dx[i+1]+k]; x[k+1:dx[i+1]]; x[dx[i+1]+j]], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00));
                    X1 = monomials([x[1:lv-dx[i+1]+k]; x[k+1:dx[i+1]]; x[dx[i+1]+j]], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1));
                    X2 = monomials([x[1:lv-dx[i+1]+k]; x[k+1:dx[i+1]]; x[dx[i+1]+j]], 0:2*1); σ2 = @variable(model, [1:1], Poly(X2));
                    X3 = monomials([x[1:lv-dx[i+1]+k]; x[k+1:dx[i+1]]; x[dx[i+1]+j]], 0:2*1); σ3 = @variable(model, [1:1], Poly(X3));
                    obj = obj - σ00[1] - σ1[1] * (x[dx[i+1]+j]*(x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])) - σ2[1] * x[dx[i+1]+j] - σ3[1] * (x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])
                    @constraints(model, begin
                        σ00 .>= 0; σ2 .>= 0; σ3 .>= 0;
                    end)
                end
            end
            for j = 1:p[i+1]-1
                for k = j+1:p[i+1]
                    for s = dx[i]:dx[i+1]-lv
                        X00 = monomials([x[s+1:s+lv]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ00 = @variable(model, [1:1], Poly(X00));
                        X1 = monomials([x[s+1:s+lv]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ1 = @variable(model, [1:1], Poly(X1));
                        obj = obj - σ00[1] + σ1[1]*((x[dx[i+1]+k] - x[dx[i+1]+j])^2 - (x[dx[i+1]+k] - x[dx[i+1]+j])*(P[i][k,:]'*x[dx[i]+1:dx[i+1]] + ξ[i][k] - P[i][j,:]'*x[dx[i]+1:dx[i+1]] - ξ[i][j]))
                        @constraints(model, begin
                            σ00 .>= 0; σ1 .>= 0
                        end)
                    end
                    for s = dx[i+1]-lv+1:dx[i+1]-1
                        X00 = monomials([x[1:lv-dx[i+1]+s]; x[s+1:dx[i+1]]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ00 = @variable(model, [1:1], Poly(X00));
                        X1 = monomials([x[1:lv-dx[i+1]+s]; x[s+1:dx[i+1]]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ1 = @variable(model, [1:1], Poly(X1));
                        obj = obj - σ00[1] + σ1[1]*((x[dx[i+1]+k] - x[dx[i+1]+j])^2 - (x[dx[i+1]+k] - x[dx[i+1]+j])*(P[i][k,:]'*x[dx[i]+1:dx[i+1]] + ξ[i][k] - P[i][j,:]'*x[dx[i]+1:dx[i+1]] - ξ[i][j]))
                        @constraints(model, begin
                            σ00 .>= 0; σ1 .>= 0
                        end)
                    end
                end
            end
        end
    elseif method == "cycle_v"
        for k = dx[1]:dx[2]-lv
            # X00 = monomials(x[k+1:k+lv], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00))
            X1 = monomials(x[k+1:k+lv], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1))
            obj = obj- σ1[1] * (x[1:dx[2]]'*Q0*x[1:dx[2]].+b0'*x[1:dx[2]].+c0)
            @constraint(model, σ1 .>= 0)
        end
        for k = dx[2]-lv+1:dx[2]-1
            # X00 = monomials(x[k+1:k+lv], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00))
            X1 = monomials(x[1:lv-dx[2]+k], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1))
            obj = obj- σ1[1] * (x[1:dx[2]]'*Q0*x[1:dx[2]].+b0'*x[1:dx[2]].+c0)
            @constraint(model, σ1 .>= 0)
        end
        for i = 1:l
            for j = 1:p[i+1]
                for k = dx[i]:dx[i+1]-lv
                    set = collect(k+1:k+lv)
                    for s = dx[i]+1:dx[i+1]
                        if !(s in set) && (length(set) != p[i])
                            set_bis = sort(vcat(set, s))
                            X00 = monomials([x[set_bis]; x[dx[i+1]+j]], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00));
                            obj += - σ00[1]
                            @constraint(model, σ00 .>= 0)
                        elseif length(set) == p[i]
                            X00 = monomials([x[set]; x[dx[i+1]+j]], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00));
                            obj += - σ00[1]
                            @constraint(model, σ00 .>= 0)
                        end
                    end
                    X1 = monomials([x[k+1:k+lv]; x[dx[i+1]+j]], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1));
                    X2 = monomials([x[k+1:k+lv]; x[dx[i+1]+j]], 0:2*1); σ2 = @variable(model, [1:1], Poly(X2));
                    X3 = monomials([x[k+1:k+lv]; x[dx[i+1]+j]], 0:2*1); σ3 = @variable(model, [1:1], Poly(X3));
                    obj = obj - σ1[1] * (x[dx[i+1]+j]*(x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])) - σ2[1] * x[dx[i+1]+j] - σ3[1] * (x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])
                    @constraints(model, begin
                        σ2 .>= 0; σ3 .>= 0;
                    end)
                end
                for k = dx[i+1]-lv+1:dx[i+1]-1
                    set = vcat(collect(1:lv-dx[i+1]+k), collect(k+1:dx[i+1]))
                    for s = dx[i]+1:dx[i+1]
                        if !(s in set) && (length(set) != p[i])
                            set_bis = sort(vcat(set, s))
                            X00 = monomials([x[set_bis]; x[dx[i+1]+j]], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00));
                            obj += - σ00[1]
                            @constraint(model, σ00 .>= 0)
                        elseif length(set) == p[i]
                            X00 = monomials([x[set]; x[dx[i+1]+j]], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00));
                            obj += - σ00[1]
                            @constraint(model, σ00 .>= 0)
                        end
                    end
                    X1 = monomials([x[1:lv-dx[i+1]+k]; x[k+1:dx[i+1]]; x[dx[i+1]+j]], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1));
                    X2 = monomials([x[1:lv-dx[i+1]+k]; x[k+1:dx[i+1]]; x[dx[i+1]+j]], 0:2*1); σ2 = @variable(model, [1:1], Poly(X2));
                    X3 = monomials([x[1:lv-dx[i+1]+k]; x[k+1:dx[i+1]]; x[dx[i+1]+j]], 0:2*1); σ3 = @variable(model, [1:1], Poly(X3));
                    obj = obj - σ00[1] - σ1[1] * (x[dx[i+1]+j]*(x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])) - σ2[1] * x[dx[i+1]+j] - σ3[1] * (x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])
                    @constraints(model, begin
                        σ2 .>= 0; σ3 .>= 0;
                    end)
                end
            end
            for j = 1:p[i+1]-1
                for k = j+1:p[i+1]
                    for s = dx[i]:dx[i+1]-lv
                        set = collect(s+1:s+lv)
                        for t = dx[i]+1:dx[i+1]
                            if !(t in set) && (length(set) != p[i])
                                set_bis = sort(vcat(set, t))
                                X00 = monomials([x[set_bis]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ00 = @variable(model, [1:1], Poly(X00));
                                obj += - σ00[1]
                                @constraint(model, σ00 .>= 0)
                            elseif length(set) == p[i]
                                X00 = monomials([x[set]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ00 = @variable(model, [1:1], Poly(X00));
                                obj += - σ00[1]
                                @constraint(model, σ00 .>= 0)
                            end
                        end
                        X1 = monomials([x[s+1:s+lv]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ1 = @variable(model, [1:1], Poly(X1));
                        obj = obj - σ00[1] + σ1[1]*((x[dx[i+1]+k] - x[dx[i+1]+j])^2 - (x[dx[i+1]+k] - x[dx[i+1]+j])*(P[i][k,:]'*x[dx[i]+1:dx[i+1]] + ξ[i][k] - P[i][j,:]'*x[dx[i]+1:dx[i+1]] - ξ[i][j]))
                        @constraint(model, σ1 .>= 0)
                    end
                    for s = dx[i+1]-lv+1:dx[i+1]-1
                        set = vcat(collect(1:lv-dx[i+1]+s), collect(s+1:dx[i+1]))
                        for t = dx[i]+1:dx[i+1]
                            if !(t in set) && (length(set) != p[i])
                                set_bis = sort(vcat(set, t))
                                X00 = monomials([x[set_bis]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ00 = @variable(model, [1:1], Poly(X00));
                                obj += - σ00[1]
                                @constraint(model, σ00 .>= 0)
                            elseif length(set) == p[i]
                                X00 = monomials([x[set]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ00 = @variable(model, [1:1], Poly(X00));
                                obj += - σ00[1]
                                @constraint(model, σ00 .>= 0)
                            end
                        end
                        X1 = monomials([x[1:lv-dx[i+1]+s]; x[s+1:dx[i+1]]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ1 = @variable(model, [1:1], Poly(X1));
                        obj = obj - σ00[1] + σ1[1]*((x[dx[i+1]+k] - x[dx[i+1]+j])^2 - (x[dx[i+1]+k] - x[dx[i+1]+j])*(P[i][k,:]'*x[dx[i]+1:dx[i+1]] + ξ[i][k] - P[i][j,:]'*x[dx[i]+1:dx[i+1]] - ξ[i][j]))
                        @constraint(model, σ1 .>= 0)
                    end
                end
            end
        end
    elseif method == "random"
        set = sort(sample(1:p[1], lv, replace=false))
        #set = sort(randperm(p[1])[1:lv])
        X00 = monomials(x[set], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00))
        X1 = monomials(x[set], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1))
        obj = obj- σ1[1] * (x[1:dx[2]]'*Q0*x[1:dx[2]].+b0'*x[1:dx[2]].+c0)
        @constraint(model, σ1 .>= 0)
        for i = 1:l
            for j = 1:p[i+1]
                set = sort(sample(1:p[i], lv, replace=false))
                X00 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00));
                X1 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1));
                X2 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*1); σ2 = @variable(model, [1:1], Poly(X2));
                X3 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*1); σ3 = @variable(model, [1:1], Poly(X3));
                obj = obj - σ00[1] - σ1[1] * (x[dx[i+1]+j]*(x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])) - σ2[1] * x[dx[i+1]+j] - σ3[1] * (x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])
                @constraints(model, begin
                    σ00 .>= 0; σ2 .>= 0; σ3 .>= 0;
                end)
            end
            for j = 1:p[i+1]-1
                for k = j+1:p[i+1]
                    set = sort(sample(1:p[i], lv, replace=false))
                    X00 = monomials([x[dx[i].+set]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ00 = @variable(model, [1:1], Poly(X00));
                    X1 = monomials([x[dx[i].+set]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ1 = @variable(model, [1:1], Poly(X1));
                    obj = obj - σ00[1] + σ1[1]*((x[dx[i+1]+k] - x[dx[i+1]+j])^2 - (x[dx[i+1]+k] - x[dx[i+1]+j])*(P[i][k,:]'*x[dx[i]+1:dx[i+1]] + ξ[i][k] - P[i][j,:]'*x[dx[i]+1:dx[i+1]] - ξ[i][j]))
                    @constraints(model, begin
                        σ00 .>= 0; σ1 .>= 0
                    end)
                end
            end
        end
    elseif method == "random_v"
        set = sort(sample(1:p[1], lv, replace=false))
        #set = sort(randperm(p[1])[1:lv])
        X00 = monomials(x[set], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00))
        X1 = monomials(x[set], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1))
        obj += - σ1[1] * (x[1:dx[2]]'*Q0*x[1:dx[2]].+b0'*x[1:dx[2]].+c0)
        @constraint(model, σ1 .>= 0)
        for i = 1:l
            for j = 1:p[i+1]
                set = sort(sample(1:p[i], lv, replace=false))
                for s = 1:p[i]
                    if !(s in set) && (length(set) != p[i])
                        set_bis = sort(vcat(set, s))
                        X00 = monomials([x[dx[i].+set_bis]; x[dx[i+1]+j]], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00));
                        obj += - σ00[1]
                        @constraint(model, σ00 .>= 0)
                    elseif length(set) == p[i]
                        X00 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00));
                        obj += - σ00[1]
                        @constraint(model, σ00 .>= 0)
                    end
                end
                X1 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1));
                X2 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*1); σ2 = @variable(model, [1:1], Poly(X2));
                X3 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*1); σ3 = @variable(model, [1:1], Poly(X3));
                obj += - σ1[1] * (x[dx[i+1]+j]*(x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])) - σ2[1] * x[dx[i+1]+j] - σ3[1] * (x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])
                @constraints(model, begin
                    σ2 .>= 0; σ3 .>= 0;
                end)
            end
            for j = 1:p[i+1]-1
                for k = j+1:p[i+1]
                    set = sort(sample(1:p[i], lv, replace=false))
                    for s = 1:p[i]
                        if !(s in set) && (length(set) != p[i])
                            set_bis = sort(vcat(set, s))
                            X00 = monomials([x[dx[i].+set_bis]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ00 = @variable(model, [1:1], Poly(X00));
                            obj += - σ00[1]
                            @constraint(model, σ00 .>= 0)
                        elseif length(set) != p[i]
                            X00 = monomials([x[dx[i].+set]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ00 = @variable(model, [1:1], Poly(X00));
                            obj += - σ00[1]
                            @constraint(model, σ00 .>= 0)
                        end
                    end
                    X1 = monomials([x[dx[i].+set]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ1 = @variable(model, [1:1], Poly(X1));
                    obj += σ1[1]*((x[dx[i+1]+k] - x[dx[i+1]+j])^2 - (x[dx[i+1]+k] - x[dx[i+1]+j])*(P[i][k,:]'*x[dx[i]+1:dx[i+1]] + ξ[i][k] - P[i][j,:]'*x[dx[i]+1:dx[i+1]] - ξ[i][j]))
                    @constraint(model, σ1 .>= 0)
                end
            end
        end
    elseif method == "order"
        set = collect(1:lv)
        #set = sort(randperm(p[1])[1:lv])
        X00 = monomials(x[set], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00))
        X1 = monomials(x[set], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1))
        obj = obj- σ1[1] * (x[1:dx[2]]'*Q0*x[1:dx[2]].+b0'*x[1:dx[2]].+c0)
        @constraint(model, σ1 .>= 0)
        for i = 1:l
            for j = 1:p[i+1]
                set = sort(sample(1:p[i], lv, replace=false))
                X00 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00));
                X1 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1));
                X2 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*1); σ2 = @variable(model, [1:1], Poly(X2));
                X3 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*1); σ3 = @variable(model, [1:1], Poly(X3));
                obj = obj - σ00[1] - σ1[1] * (x[dx[i+1]+j]*(x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])) - σ2[1] * x[dx[i+1]+j] - σ3[1] * (x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])
                @constraints(model, begin
                    σ00 .>= 0; σ2 .>= 0; σ3 .>= 0;
                end)
            end
            for j = 1:p[i+1]-1
                for k = j+1:p[i+1]
                    set = sort(sample(1:p[i], lv, replace=false))
                    X00 = monomials([x[dx[i].+set]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ00 = @variable(model, [1:1], Poly(X00));
                    X1 = monomials([x[dx[i].+set]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ1 = @variable(model, [1:1], Poly(X1));
                    obj = obj - σ00[1] + σ1[1]*((x[dx[i+1]+k] - x[dx[i+1]+j])^2 - (x[dx[i+1]+k] - x[dx[i+1]+j])*(P[i][k,:]'*x[dx[i]+1:dx[i+1]] + ξ[i][k] - P[i][j,:]'*x[dx[i]+1:dx[i+1]] - ξ[i][j]))
                    @constraints(model, begin
                        σ00 .>= 0; σ1 .>= 0
                    end)
                end
            end
        end
    elseif method == "order_v"
        set = collect(1:lv)
        #set = sort(randperm(p[1])[1:lv])
        X00 = monomials(x[set], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00))
        X1 = monomials(x[set], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1))
        obj += - σ1[1] * (x[1:dx[2]]'*Q0*x[1:dx[2]].+b0'*x[1:dx[2]].+c0)
        @constraint(model, σ1 .>= 0)
        for i = 1:l
            for j = 1:p[i+1]
                for s = 1:p[i]
                    if !(s in set) && (length(set) != p[i])
                        set_bis = sort(vcat(set, s))
                        X00 = monomials([x[dx[i].+set_bis]; x[dx[i+1]+j]], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00));
                        obj += - σ00[1]
                        @constraint(model, σ00 .>= 0)
                    elseif length(set) == p[i]
                        X00 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*2); σ00 = @variable(model, [1:1], Poly(X00));
                        obj += - σ00[1]
                        @constraint(model, σ00 .>= 0)
                    end
                end
                X1 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*1); σ1 = @variable(model, [1:1], Poly(X1));
                X2 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*1); σ2 = @variable(model, [1:1], Poly(X2));
                X3 = monomials([x[dx[i].+set]; x[dx[i+1]+j]], 0:2*1); σ3 = @variable(model, [1:1], Poly(X3));
                obj += - σ1[1] * (x[dx[i+1]+j]*(x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])) - σ2[1] * x[dx[i+1]+j] - σ3[1] * (x[dx[i+1]+j]-P[i][j,:]'*x[dx[i]+1:dx[i+1]]-ξ[i][j])
                @constraints(model, begin
                    σ2 .>= 0; σ3 .>= 0;
                end)
            end
            for j = 1:p[i+1]-1
                for k = j+1:p[i+1]
                    for s = 1:p[i]
                        if !(s in set) && (length(set) != p[i])
                            set_bis = sort(vcat(set, s))
                            X00 = monomials([x[dx[i].+set_bis]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ00 = @variable(model, [1:1], Poly(X00));
                            obj += - σ00[1]
                            @constraint(model, σ00 .>= 0)
                        elseif length(set) != p[i]
                            X00 = monomials([x[dx[i].+set]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ00 = @variable(model, [1:1], Poly(X00));
                            obj += - σ00[1]
                            @constraint(model, σ00 .>= 0)
                        end
                    end
                    X1 = monomials([x[dx[i].+set]; x[dx[i+1]+j]; x[dx[i+1]+k]], 0:2*0); σ1 = @variable(model, [1:1], Poly(X1));
                    obj += σ1[1]*((x[dx[i+1]+k] - x[dx[i+1]+j])^2 - (x[dx[i+1]+k] - x[dx[i+1]+j])*(P[i][k,:]'*x[dx[i]+1:dx[i+1]] + ξ[i][k] - P[i][j,:]'*x[dx[i]+1:dx[i+1]] - ξ[i][j]))
                    @constraint(model, σ1 .>= 0)
                end
            end
        end
    end
    @constraint(model, obj .== 0)
    # solve
    MOI.set(model, MOI.Silent(), true);
    optimize!(model)
    @printf("Optimal value: %.2f, termination status: %s, primal status: %s, dual status: %s, solving time: %.2f.\n", value.(lower_bound), termination_status(model), primal_status(model), dual_status(model), solve_time(model))
    # println(value.(Q)); println(value.(b)); println(value.(c))
    # return value.(-Q*Q), value.(-2 .* Q*b), value.(1 .- b'*b)
    return value.(Q), value.(b), value.(c), value.(lower_bound)
end

# compute the outer approximation ellipsoid
function OuterApproximationMorari(Q0, b0, c0, P, ξ)
    l = length(P);
    p = Int.(zeros(l+1,1)); p[1] = size(P[1],2); d = 0;
    for i = 1:l
        p[i+1] = size(P[i],1)
        d = d + p[i]
    end
    model = Model(with_optimizer(Mosek.Optimizer))
    # objective variable
    Q = @variable(model, [1:p[l+1], 1:p[l+1]], Symmetric);
    b = @variable(model, [1:p[l+1], 1:1]);
    c = @variable(model, [1:1, 1:1]);
    μ = @variable(model, [1:1, 1:1]);
    ν = @variable(model, [1:d-p[1]+p[l+1], 1:1]);
    η = @variable(model, [1:d-p[1]+p[l+1], 1:1]);
    λ = @variable(model, [1:d, 1:1]);
    Λ = @variable(model, [1:d, 1:d]);
    # matrix components
    M1 = μ.*[Q0 1/2*b0; 1/2*b0' c0];
    B_M1 = [Matrix(I(p[1])) zeros(p[1],d-p[1]+p[l+1]+1); zeros(1,d+p[l+1]) 1];
    M1 = B_M1'*M1*B_M1;
    #
    E = Matrix(I(d-p[1]+p[l+1])); T = zeros(d-p[1]+p[l+1],d-p[1]+p[l+1]);
    for i = 1:d-p[1]+p[l+1]
        T = T + λ[i].*E[i,:]*E[i,:]'
    end
    for i = 1:d-p[1]+p[l+1]-1
        for j = i+1:d-p[1]+p[l+1]
            T = T + Λ[i,j].*(E[i,:] - E[j,:])*(E[i,:] - E[j,:])'
        end
    end
    M2 = [zeros(d-p[1]+p[l+1],d-p[1]+p[l+1]) T -ν; T -2*T ν+η; -ν' ν'+η' 0];
    B_M2 = [P[1] zeros(p[2],sum(p[2:end])) ξ[1]];
    for i = 1:l-1
        B_M2 = vcat(B_M2, [zeros(p[i+2],sum(p[1:i])) P[i+1] zeros(p[i+2],sum(p[i+2:end])) ξ[i+1]])
    end
    B_M2 = vcat(B_M2, [zeros(d-p[1]+p[l+1],p[1]) Matrix(I(d-p[1]+p[l+1])) zeros(d-p[1]+p[l+1],1); zeros(1,d+p[l+1]) 1])
    M2 = B_M2'*M2*B_M2
    #
    # S = [-Q -1/2*b; -1/2*b' -c[1]];
    # B_S = [zeros(p[l+1],d) P[l] zeros(p[l+1],1); zeros(1,d+p[l+1]) 1];
    # S = B_S'*S*B_S
    e = [zeros(d+p[l+1],1); 1];
    M = [[M1+M2-e*e' [zeros(d,p[l+1]); Matrix(I(p[l+1]))*Q; b']]; [zeros(p[l+1],d) Q*Matrix(I(p[l+1])) b -Matrix(I(p[l+1]))]];
    @variable(model, lower_bound)
    @objective(model, Max, lower_bound)
    vect = lower_bound
    for i = 1:p[l+1]
        vect = hcat(vect, -Q[i,i:p[l+1]]')
    end
    @constraints(model, begin
        μ .>= 0
        ν .>= 0
        η .>= 0
        Λ .>= 0
        vect[:] in MOI.RootDetConeTriangle(p[l+1])
        -M in PSDCone()
    end)
    # -M in PSDCone()
    # -M1-M2-M3-M4-M5 in PSDCone()
    # solve
    MOI.set(model, MOI.Silent(), true);
    optimize!(model)
    @printf("Optimal value: %.2f, termination status: %s, primal status: %s, dual status: %s, solving time: %.2f.\n", value.(lower_bound.^2), termination_status(model), primal_status(model), dual_status(model), solve_time(model))
    # println(value.(Q)); println(value.(b)); println(value.(c))
    return value.(-Q*Q), value.(-2 .* Q*b), value.(1 .- b'*b), value.(lower_bound.^2)
    # return value.(Q), value.(b), value.(c)
end

function OuterApproximationFazlyab(Q0, b0, c0, P, ξ)
    p1 = size(P[1], 2); p2 = size(P[1], 1)
    model = Model(with_optimizer(Mosek.Optimizer))
    # objective variable
    Q = @variable(model, [1:p2, 1:p2], Symmetric);
    b = @variable(model, [1:p2, 1:1]);
    c = @variable(model, [1:1, 1:1]);
    λ1 = @variable(model, [1:1, 1:1]);
    λ2 = @variable(model, [1:p2, 1:1]);
    λ3 = @variable(model, [1:p2, 1:1]);
    λ4 = @variable(model, [1:p2, 1:1]);
    Λ = @variable(model, [1:p2, 1:p2]);
    # matrix components
    M1 = [zeros(p1,p1+p2+1); zeros(p2,p1) -Q -1/2*b; zeros(1,p1) -1/2*b' -c[1]];
    M2 = λ1.*[Q0 zeros(p1,p2) 1/2*b0; zeros(p2,p1+p2+1); 1/2*b0' zeros(1,p2) c0];
    M3 = [zeros(p1,p1) -1/2*P[1]'*diagm(λ2[:]) zeros(p1,1); -1/2*diagm(λ2[:])*P[1] diagm(λ2[:]) -1/2*(λ2 .* ξ[1]); zeros(1,p1) -1/2*(λ2' .* ξ[1]') 0];
    M4 = [zeros(p1,p1+p2+1); zeros(p2,p1+p2) 1/2*λ3; zeros(1,p1) 1/2*λ3' 0];
    M5 = [zeros(p1,p1+p2) -1/2*P[1]'*λ4; zeros(p2,p1+p2) 1/2*λ4; -1/2*λ4'*P[1] 1/2*λ4' -λ4'*ξ[1]];
    E = Matrix(I(p2)); T = zeros(p2,p2);
    for i = 1:p2-1
        for j = i+1:p2
            T = T + Λ[i,j].*(E[i,:] - E[j,:])*(E[i,:] - E[j,:])'
        end
    end
    M = M1+M2+M3+M4+M5+[zeros(p1,p1) 1/2*P[1]'*T zeros(p1,1); 1/2*T*P[1] -T 1/2*T*ξ[1]; zeros(1,p1) 1/2*ξ[1]'*T 0]
    # M = [[M2+M3+M4+M5+[zeros(p1,p1) 1/2*P[1]'*T zeros(p1,1); 1/2*T*P[1] -T 1/2*T*ξ[1]; zeros(1,p1) 1/2*ξ[1]'*T 0]-e*e' [zeros(p1,p2); I(p2)*Q; b']]; [zeros(p2,p1) Q*I(p2) b -Matrix(I(p2))]];
    @variable(model, lower_bound)
    @objective(model, Max, lower_bound)
    Mat = [-Q b/2; b'/2 1-c[1]];
    vect = lower_bound
    for i = 1:p2
        vect = hcat(vect, -Q[i,i:p2]')
    end
    @constraints(model, begin
        λ1 .>= 0
        λ3 .>= 0
        λ4 .>= 0
        Λ .>= 0
        Mat in PSDCone()
        vect[:] in MOI.RootDetConeTriangle(p2)
        -M in PSDCone()
    end)
    # -M in PSDCone()
    # -M1-M2-M3-M4-M5 in PSDCone()
    # solve
    MOI.set(model, MOI.Silent(), true);
    optimize!(model)
    @printf("Optimal value: %.2f, termination status: %s, primal status: %s, dual status: %s.\n", value.(lower_bound.^2), termination_status(model), primal_status(model), dual_status(model))
    # println(value.(Q)); println(value.(b)); println(value.(c))
    # return value.(-Q*Q), value.(-2 .* Q*b), value.(1 .- b'*b)
    return value.(Q), value.(b), value.(c), value.(lower_bound)
end

# plot the ellipsoid
function OuterApproximationPlot(Q0, b0, c0, ord, P, ξ, method, k)
    θ = 0:0.1:2*π+0.1
    D = svd(P[1]); Q = D.U*diagm(D.S)^(-1)*D.V'*Q0*D.V*diagm(D.S)^(-1)*D.U'; b = D.U*diagm(D.S)^(-1)*D.V'*b0; c = c0
    F = eigen(Q); T = F.vectors; Γ = diagm(sqrt.(-F.values));
    v = -1/2 * Q^(-1) * b; m = c .- v'*Q*v
    x11 = ((sqrt(m).*T*Γ^(-1)*[cos.(θ)'; sin.(θ)'].+v).+ξ[1])[1,:] .* (((sqrt(m).*T*Γ^(-1)*[cos.(θ)'; sin.(θ)'].+v).+ξ[1])[1,:] .>= 0);
    x21 = ((sqrt(m).*T*Γ^(-1)*[cos.(θ)'; sin.(θ)'].+v).+ξ[1])[2,:] .* (((sqrt(m).*T*Γ^(-1)*[cos.(θ)'; sin.(θ)'].+v).+ξ[1])[2,:] .>= 0);
    for i = 2:k
        X = [x11'; x21']
        x11 = (P[i]*X.+ξ[i])[1,:] .* ((P[i]*X.+ξ[i])[1,:] .>= 0);
        x21 = (P[i]*X.+ξ[i])[2,:] .* ((P[i]*X.+ξ[i])[2,:] .>= 0);
    end
    Q = Q0; b = b0; c = c0;
    for i = 1:k
        Q, b, c = OuterApproximation(Q, b, c, ord, P[i], ξ[i], method)
    end
    # @printf("Q = "); println(Q); @printf("b = "); println(b);  @printf("c = %f\n", c)
    # @printf("Determinant: %f\n", det(-Q))
    F = eigen(Q); T = F.vectors; Γ = diagm(sqrt.(-F.values));
    v = -1/2 * Q^(-1) * b; m = c - v'*Q*v
    x12 = (sqrt(m).*T*Γ^(-1)*[cos.(θ)'; sin.(θ)'].+v)[1,:];
    x22 = (sqrt(m).*T*Γ^(-1)*[cos.(θ)'; sin.(θ)'].+v)[2,:];
    x1 = [x11 x12]; x2 = [x21 x22];
    p = plot(x1, x2, title = @sprintf("Ord %d, %s", ord, uppercasefirst(method)), label = ["Image" "OutApprox"], legend=false)#:outertopright)
    return p
end

# plot the ellipsoid by sampling
function OuterApproximationPlotSampling(Q0, b0, c0, ord, P, ξ, method, k; lv=[], morari=[], meth="cycle")
    num = 100000; n = size(Q0, 1); x11 = zeros(num,1); x21 = zeros(num,1);
    Q = Q0; b = b0; c = c0;
    F = eigen(Q); T = F.vectors; Γ = diagm(sqrt.(-F.values));
    v = -1/2 * Q^(-1) * b; m = c .- v'*Q*v
    for i = 1:num
        x = randn(n,1); x = x./norm(x, 2); x = sqrt(m).*T*Γ^(-1)*x.+v
        for j = 1:k
            x = (P[j]*x+ξ[j]) .* (P[j]*x+ξ[j] .>= 0)
        end
        x11[i] = x[1]; x21[i] = x[2]
    end
    if method == "Morari"
        Q, b, c, mo = OuterApproximationMorari(Q, b, c, P[1:k], ξ[1:k])
    elseif method == "Fazlyab"
        Q, b, c = OuterApproximationFazlyab(Q, b, c, P[1:k], ξ[1:k])
    elseif method == "sublevel"
        Q, b, c, sub = OuterApproximationSublevel(Q, b, c, P[1:k], ξ[1:k], lv, method=meth)
    else
        for i = 1:k
            Q, b, c = OuterApproximation(Q, b, c, ord, P[i], ξ[i], method)
        end
    end
    # @printf("Q = "); println(Q); @printf("b = "); println(b);  @printf("c = %f\n", c)
    # @printf("Determinant: %f\n", det(-Q))
    F = eigen(Q); T = F.vectors; Γ = diagm(sqrt.(-F.values));
    v = -1/2 * Q^(-1) * b; m = c .- v'*Q*v
    θ = 0:0.1:2*π+0.1;
    x12 = (sqrt(m).*T*Γ^(-1)*[cos.(θ)'; sin.(θ)'].+v)[1,:];
    x22 = (sqrt(m).*T*Γ^(-1)*[cos.(θ)'; sin.(θ)'].+v)[2,:];
    # x1 = [x11 x12]; x2 = [x21 x22];
    p = plot(x11, x21, seriestype = :scatter, markersize = 1)
    if method in ["Morari" "Fazlyab"]
        p = plot!(x12, x22, title = method, label = ["Image" "OutApprox"], legend=false)#:outertopright)
    elseif method == "sublevel"
        p = plot!(x12, x22, title = @sprintf("d=%d, lv=%d, r=%.2f%%", ord, lv, (sub[1]-morari[1])/morari[1] * 100), label = ["Image" "OutApprox"], legend=false)
    else
        p = plot!(x12, x22, title = @sprintf("d=%d, %s", ord, uppercasefirst(method)), label = ["Image" "OutApprox"], legend=false)#:outertopright)
    end
    if method == "Morari"
        return p, mo
    else
        return p
    end
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

# n = 2; Q01 = Matrix(-I(n)); b01 = zeros(n,1); c01 = 1;
# L1 = 1; P1 = Array{Any}(undef, L1); ξ1 = Array{Any}(undef, L1);
# for i = 1:L1-1
#     P1[i] = rand(n,n); ξ1[i] = rand(n,1);
#     for j = 1:n
#         P1[i][j,j] = 1
#     end
# end
# P1[L1] = rand(2,n); ξ1[L1] = rand(2,1); P1[L1][1,1] = 1; P1[L1][2,2] = 1;
# # ξ1[L1] = zeros(2,1)
#
# # plot original graph and image under ReLU
# p11 = OuterApproximationPlotSampling(Q01, b01, c01, 1, P1, ξ1, "tong", L1);
# p12 = OuterApproximationPlotSampling(Q01, b01, c01, 2, P1, ξ1, "tong", L1);
# p21 = OuterApproximationPlotSampling(Q01, b01, c01, 1, P1, ξ1, "split", L1);
# p22 = OuterApproximationPlotSampling(Q01, b01, c01, 2, P1, ξ1, "split", L1);
# p31, m = OuterApproximationPlotSampling(Q01, b01, c01, 1, P1, ξ1, "Morari", L1);
# p32 = OuterApproximationPlotSampling(Q01, b01, c01, 2, P1, ξ1, "sublevel", L1, lv=2, morari=m, meth="random");
#
# plot(p11, p21, p31, p12, p22, p32, layout = grid(2,3), fmt = :png)





# n = 10; Q01 = Matrix(-I(n)); b01 = zeros(n,1); c01 = 1;
# L1 = 2; P1 = Array{Any}(undef, L1); ξ1 = Array{Any}(undef, L1);
# for i = 1:L1-1
#     P1[i] = randn(n,n); ξ1[i] = randn(n,1);
#     for j = 1:n
#         P1[i][j,j] = 1
#     end
# end
# P1[L1] = rand(2,n); ξ1[L1] = rand(2,1); P1[L1][1,1] = 1; P1[L1][2,2] = 1;
#
# p11 = plot();
# p21, m = OuterApproximationPlotSampling(Q01, b01, c01, 1, P1, ξ1, "Morari", L1);
# p31 = plot();
# p12 = plot();
# p22 = OuterApproximationPlotSampling(Q01, b01, c01, 2, P1, ξ1, "sublevel", L1, lv=2, morari=m);
# p32 = plot();
# p13 = OuterApproximationPlotSampling(Q01, b01, c01, 2, P1, ξ1, "sublevel", L1, lv=3, morari=m);
# p23 = OuterApproximationPlotSampling(Q01, b01, c01, 2, P1, ξ1, "sublevel", L1, lv=5, morari=m);
# p33 = OuterApproximationPlotSampling(Q01, b01, c01, 2, P1, ξ1, "sublevel", L1, lv=7, morari=m);
# p23 = plot(); p33 = plot();
#
# plot(p11, p21, p31, p12, p22, p32, p13, p23, p33, layout = grid(3,3), fmt = :png)





n = 20; Q10 = Matrix(-I(n)); b10 = zeros(n,1); c10 = 1;
L10 = 2; P10 = Array{Any}(undef, L10); ξ10 = Array{Any}(undef, L10);
for i = 1:L10-1
    P10[i] = randn(n,n)*0.5; ξ10[i] = randn(n,1);
    for j = 1:n
        P10[i][j,j] = 1
    end
end
P10[L10] = rand(2,n)*0.5; ξ10[L10] = rand(2,1); P10[L10][1,1] = 1; P10[L10][2,2] = 1;

# plot original graph and image under ReLU
p11, m = OuterApproximationPlotSampling(Q10, b10, c10, 1, P10, ξ10, "Morari", L10);
# p21 = OuterApproximationPlotSampling(Q10, b10, c10, 1, P10, ξ10, "jean", L10);
p21 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=2, morari=m, meth="cycle_v");
p31 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=2, morari=m, meth="order_v");
p12 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=3, morari=m, meth="order_v");
# p22 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "jean", L10);
p22 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=4, morari=m, meth="order_v");
p32 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=5, morari=m, meth="order_v");
p13 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=6, morari=m, meth="order_v");
p23 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=7, morari=m, meth="order_v");
p33 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=8, morari=m, meth="order_v");

plot(p11, p21, p31, p12, p22, p32, p13, p23, p33, layout = Plots.grid(3,3), fmt = :png)


# A = zeros(9,9);
# for i = 1:3
#     for j = i:6
#         A[i,j] = 1; A[j,i] = 1;
#     end
# end
# for i = 4:9
#     for j = 7:9
#         A[i,j] = 1; A[j,i] = 1;
#     end
# end
# cliques, NumElem, A_new, g, g_new = ChordalExtension(A, turn="on")
