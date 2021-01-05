function solve_moment_lip(A, b, c, x00, e, options; cliques = [], sizes = [])
    if isequal(options["range"], "global")
        ϵ = 10;
    else
        ϵ = e;
    end
    d0 = size(A[1], 2); d = size(A[1], 1); n = 2*d0+d; typ = "max";
    lx = x00 .- ϵ*ones(d0,1); ux = x00 .+ ϵ*ones(d0,1);
    lz = (A[1] .* (A[1] .> 0)) * lx .+ (A[1] .* (A[1] .< 0)) * ux .+ b[1];
    uz = (A[1] .* (A[1] .> 0)) * ux .+ (A[1] .* (A[1] .< 0)) * lx .+ b[1];
    ly = 1*(lz .> 0); uy = 1*(uz .> 0);
    E = sparse([I(n); I(n)]);
    obj = [zeros(1+d0,n+1); [zeros(d,1+d0+d) 1/2*diagm(c[:])*A[1]; zeros(d0,1+d0) 1/2*A[1]'*diagm(c[:]) zeros(d0,d0)]];
    if options["level"] == 0 && options["clique"] == "on"# sparse Shor's relaxation
        MomConst = Array{Any}(undef, 0); NumMom = 0
        for i = 1:length(sizes)
            NumMom += 1
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[cliques[i]', :]);
            MomConst[NumMom]["ord"] = 1;
        end
        LocConst = Array{Any}(undef, 0); NumLoc = 0
        for i = 1:d0
            # t^2<=1
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[d0+d+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = 0;
            # (x-lx)(x-ux)<=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([lx[i]*ux[i]*100 zeros(1,n); -(lx[i]+ux[i])*10 E[i,:]'; 1*1 2*E[i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "<=";
            LocConst[NumLoc]["ord"] = 0;
        end
        for i = 1:d
            #y*(y-1)==0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[d0+i,:]'; -1 2*E[d0+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "==";
            LocConst[NumLoc]["ord"] = 0;
            #(y-ly)(y-uy)<=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([ly[i]*uy[i] zeros(1,n); -(ly[i]+uy[i]) E[d0+i,:]'; 1 2*E[d0+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "<=";
            LocConst[NumLoc]["ord"] = 0;
            #(y-1/2)*(Ax+b)>=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([-1/2*b[1][i]*10 zeros(1,n); b[1][i]*10 E[d0+i,:]'; [-1/2*A[1][i,:] E[1:d0, :]]; [A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1: 2*d0+d]]]);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = 0;
        end
    elseif options["level"] == 0 && options["clique"] == "off"# dense Shor's relaxation
        MomConst = Array{Any}(undef, 0); NumMom = 0
        NumMom += 1
        MomConst = vcat(MomConst, Array{Any}(undef, 1));
        MomConst[NumMom] = Dict();
        MomConst[NumMom]["basis"] = sparse(I(n));
        MomConst[NumMom]["ord"] = 1;
        LocConst = Array{Any}(undef, 0); NumLoc = 0
        for i = 1:d0
            # t^2<=1
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[d0+d+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = 0;
            # (x-lx)(x-ux)<=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([lx[i]*ux[i]*100 zeros(1,n); -(lx[i]+ux[i])*10 E[i,:]'; 1*1 2*E[i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "<=";
            LocConst[NumLoc]["ord"] = 0;
        end
        for i = 1:d
            #y*(y-1)==0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[d0+i,:]'; -1 2*E[d0+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "==";
            LocConst[NumLoc]["ord"] = 0;
            #(y-ly)(y-uy)<=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([ly[i]*uy[i] zeros(1,n); -(ly[i]+uy[i]) E[d0+i,:]'; 1 2*E[d0+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "<=";
            LocConst[NumLoc]["ord"] = 0;
            #(y-1/2)*(Ax+b)>=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([-1/2*b[1][i]*10 zeros(1,n); b[1][i]*10 E[d0+i,:]'; [-1/2*A[1][i,:] E[1:d0, :]]; [A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1: 2*d0+d]]]);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = 0;
        end
    elseif options["level"] > 0 && options["clique"] == "on"# sparse sublevel relaxation
        MomConst = Array{Any}(undef, 0); NumMom = 0
        for i = 1:length(sizes)
            NumMom += 1
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[cliques[i]', :]);
            MomConst[NumMom]["ord"] = 1;
        end
        LocConst = Array{Any}(undef, 0); NumLoc = 0
        for i = 1:d0
            if options["level"] < d+1
                set = filter(x->x≠(d0+d+i), cliques[d0-i+1]');
                set = [set; set];
                subset = vcat(d0+d+i, set[1:1+options["level"]-2])
            else
                subset = vcat(d0+d+i, collect(d0+1:d0+d))
            end
            NumMom += 1
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[sort(subset), :]);
            MomConst[NumMom]["ord"] = options["ord"];
            # t^2<=1
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[d0+d+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = options["ord"] - 1;
            set1 = filter(x->x≠i, collect(1:d0));
            set1 = [set1; set1];
            set2 = [collect(d0+1:d0+d); collect(d0+1:d0+d)];
            subset = vcat(i, set1[i:i+Int(options["level"]/2)-2], set2[1:1+Int(options["level"]/2)-1])
            NumMom += 1
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[sort(subset), :]);
            MomConst[NumMom]["ord"] = options["ord"];
            # (x-lx)(x-ux)<=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([lx[i]*ux[i] zeros(1,n); -(lx[i]+ux[i]) E[i,:]'; 1 2*E[i,:]']); #(x-lx)(x-ux)<=0
            LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
            LocConst[NumLoc]["typ"] = "<=";
            LocConst[NumLoc]["ord"] = options["ord"] - 1;
        end
        for i = 1:d
            set1 = [collect(1:d0); collect(1:d0)];
            set2 = filter(x->x≠(d0+i), collect(d0+1:d0+d));
            set2 = [set2; set2];
            subset = vcat(d0+i, set1[1:1+Int(options["level"]/2)-1], set2[i:i+Int(options["level"]/2)-2])
            NumMom += 1
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[sort(subset), :]);
            MomConst[NumMom]["ord"] = options["ord"];
            # y*(y-1)==0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[d0+i,:]'; -1 2*E[d0+i,:]']); #y*(y-1)==0
            LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
            LocConst[NumLoc]["typ"] = "==";
            LocConst[NumLoc]["ord"] = options["ord"] - 1;
            # (y-ly)(y-uy)<=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([ly[i]*uy[i] zeros(1,n); -(ly[i]+uy[i]) E[d0+i,:]'; 1 2*E[d0+i,:]']); #(y-ly)(y-uy)<=0
            LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
            LocConst[NumLoc]["typ"] = "<=";
            LocConst[NumLoc]["ord"] = options["ord"] - 1;
            # (y-1/2)*(Ax+b)>=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([-1/2*b[1][i] zeros(1,n); b[1][i] E[d0+i,:]'; [-1/2*A[1][i,:] E[1:d0, :]]; [A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1: 2*d0+d]]]);
            LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = options["ord"] - 1;
            if options["depth"] > 1
                depth = options["depth"]
                idx = [1+Int(options["level"]/2)*j for j in 1:depth-1];
                for k = 1:depth-1
                    set1 = [collect(1:d0); collect(1:d0)];
                    set2 = filter(x->x≠(d0+i), collect(d0+1:d0+d));
                    set2 = [set2; set2];
                    subset = vcat(d0+i, set1[idx[k]:idx[k]+Int(options["level"]/2)-1], set2[i+idx[k]-1:i+idx[k]-1+Int(options["level"]/2)-2])
                    set = [collect(1:d0); collect(1:d0)];
                    subset = vcat(d0+i, set[idx[k]:idx[k]+options["level"]-2])
                    NumMom += 1
                    MomConst = vcat(MomConst, Array{Any}(undef, 1));
                    MomConst[NumMom] = Dict();
                    MomConst[NumMom]["basis"] = sparse(E[sort(subset), :]);
                    MomConst[NumMom]["ord"] = options["ord"];
                    NumLoc += 1
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] = sparse([-1/2*b[1][i] zeros(1,n); b[1][i] E[d0+i,:]'; [-1/2*A[1][i,:] E[1:d0, :]]; [A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1: 2*d0+d]]]);
                    LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
                    LocConst[NumLoc]["typ"] = ">=";
                    LocConst[NumLoc]["ord"] = options["ord"] - 1;
                end
            end
        end
    else
        error("NotImplementedError")
    end
    @printf("\n%s LCEP problem: clique %s, order %d, level %d, depth %d\n", uppercasefirst(options["range"]), uppercase(options["clique"]), options["ord"], options["level"], options["depth"]);
    OptVal, running_time, stat = solve_moment(typ, obj, MomConst, LocConst, options; cliques = cliques, sizes = sizes);
    return OptVal, running_time, stat
end

function solve_moment_cert(A, b, c, x00, eps, options; cliques = [], sizes = [])
    if isequal(options["range"], "global")
        ϵ = 10;
    else
        ϵ = eps;
    end
    d0 = size(A[1], 2); d = size(A[1], 1); n = d0+d; typ = "max";
    lx = x00 .- ϵ*ones(d0,1); ux = x00 .+ ϵ*ones(d0,1);
    lz1 = (A[1] .* (A[1] .> 0)) * lx .+ (A[1] .* (A[1] .< 0)) * ux .+ b[1];
    uz1 = (A[1] .* (A[1] .> 0)) * ux .+ (A[1] .* (A[1] .< 0)) * lx .+ b[1];
    lx1 = lz1.*(lz1.>=0); ux1 = uz1.*(uz1.>=0);
    E = sparse([I(n); I(n)]);
    obj = [zeros(1,d0+1) 1/2*c'; zeros(d0,1+n); 1/2*c zeros(d,n)];
    if options["level"] == 0 && options["clique"] == "on"# sparse Shor's relaxation
        MomConst = Array{Any}(undef, 0); NumMom = 0
        for i = 1:length(sizes)
            NumMom += 1
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[cliques[i]', :]);
            MomConst[NumMom]["ord"] = 1;
        end
        LocConst = Array{Any}(undef, 0); NumLoc = 0
        for i = 1:d0
            #(x-lx)(x-ux)<=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([lx[i]*ux[i] zeros(1,n); -lx[i]-ux[i] E[i,:]'; 1 2*E[i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "<=";
            LocConst[NumLoc]["ord"] = 0;
        end
        for i = 1:d
            #(x1-lx1)(x1-ux1)<=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([lx1[i]*ux1[i] zeros(1,n); -lx1[i]-ux1[i] E[d0+i,:]'; 1 2*E[d0+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "<=";
            LocConst[NumLoc]["ord"] = 0;
            #x1>=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[d0+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = 0;
            #x1>=Ax+b
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[d0+i,:]'; -A[1][i,:] E[1:d0, :]; -b[1][i] zeros(1,n)]);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = 0;
            #x1(x1-Ax-b)=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 2*E[d0+i,:]'; [-A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1:d0+d]]; -b[1][i] E[d0+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "==";
            LocConst[NumLoc]["ord"] = 0;
        end
    elseif options["level"] == 0 && options["clique"] == "off"# dense Shor's relaxation
        MomConst = Array{Any}(undef, 0); NumMom = 0
        NumMom += 1;
        MomConst = vcat(MomConst, Array{Any}(undef, 1));
        MomConst[NumMom] = Dict();
        MomConst[NumMom]["basis"] = sparse(I(n));
        MomConst[NumMom]["ord"] = 1;
        LocConst = Array{Any}(undef, 0); NumLoc = 0
        for i = 1:d0
            #(x-lx)(x-ux)<=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([lx[i]*ux[i] zeros(1,n); -lx[i]-ux[i] E[i,:]'; 1 2*E[i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "<=";
            LocConst[NumLoc]["ord"] = 0;
        end
        for i = 1:d
            #(x1-lx1)(x1-ux1)<=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([lx1[i]*ux1[i] zeros(1,n); -lx1[i]-ux1[i] E[d0+i,:]'; 1 2*E[d0+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "<=";
            LocConst[NumLoc]["ord"] = 0;
            #x1>=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[d0+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = 0;
            #x1>=Ax+b
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[d0+i,:]'; -A[1][i,:] E[1:d0, :]; -b[1][i] zeros(1,n)]);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = 0;
            #x1(x1-Ax-b)=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 2*E[d0+i,:]'; [-A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1:d0+d]]; -b[1][i] E[d0+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "==";
            LocConst[NumLoc]["ord"] = 0;
        end
    elseif options["level"] > 0 && options["clique"] == "on" # sparse sublevel relaxation
        MomConst = Array{Any}(undef, 0); NumMom = 0
        for i = 1:length(sizes)
            NumMom += 1
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[cliques[i]', :]);
            MomConst[NumMom]["ord"] = 1;
        end
        LocConst = Array{Any}(undef, 0); NumLoc = 0
        for i = 1:d0
            if options["level"] < d0+1
                set = filter(x->x≠i, collect(1:d0));
                set = [set; set];
                subset = vcat(i, set[1:1+options["level"]-2])
            else
                subset = vcat(collect(1:d0), d0+1);
            end
            NumMom += 1
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[sort(subset), :]);
            MomConst[NumMom]["ord"] = options["ord"];
            #(x-lx)(x-ux)<=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([lx[i]*ux[i] zeros(1,n); -lx[i]-ux[i] E[i,:]'; 1 2*E[i,:]']);
            LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
            LocConst[NumLoc]["typ"] = "<=";
            LocConst[NumLoc]["ord"] = options["ord"] - 1;
        end
        for i = 1:d
            if options["level"] < d0+1
                set = collect(1:d0);
                set = [set; set];
                subset = vcat(d0+i, set[1:1+options["level"]-2])
            else
                subset = vcat(d0+i, collect(1:d0))
            end
            NumMom += 1
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[sort(subset), :]);
            MomConst[NumMom]["ord"] = options["ord"];
            #(x1-lx1)(x1-ux1)<=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([lx1[i]*ux1[i] zeros(1,n); -lx1[i]-ux1[i] E[d0+i,:]'; 1 2*E[d0+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
            LocConst[NumLoc]["typ"] = "<=";
            LocConst[NumLoc]["ord"] = options["ord"] - 1;
            #x1>=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[d0+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = options["ord"] - 1;
            #x1>=Ax+b
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[d0+i,:]'; -A[1][i,:] E[1:d0, :]; -b[1][i] zeros(1,n)]);
            LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = options["ord"] - 1;
            #x1(x1-Ax-b)=0
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 2*E[d0+i,:]'; [-A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1:d0+d]]; -b[1][i] E[d0+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
            LocConst[NumLoc]["typ"] = "==";
            LocConst[NumLoc]["ord"] = options["ord"] - 1;
        end
    else
        error("NotImplementedError")
    end
    @printf("\n%s Cert problem: clique %s, order %d, level %d, depth %d\n", uppercasefirst(options["range"]), uppercase(options["clique"]), options["ord"], options["level"], options["level"])
    OptVal, running_time, stat = solve_moment(typ, obj, MomConst, LocConst, options; cliques=cliques, sizes=sizes);
    return OptVal, running_time, stat
end

function solve_moment_maxcut(A, W, options; cliques = [], sizes = [], M = [])
    n = size(A, 1); typ = "max";
    L = sparse(diagm(((A.*W)*ones(n,1))[:,1])) - A.*W;
    E = sparse([I(n); I(n)]);
    obj = [zeros(1,n+1); [zeros(n,1) 1/4*L]];
    if options["level"] == 0 && isequal(options["clique"], "on") # sparse Shor's relaxation
        # define moment matrices
        MomConst = Array{Any}(undef, 0); NumMom = 0
        for i = 1:length(sizes)
            NumMom += 1
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[cliques[i]', :]);
            MomConst[NumMom]["ord"] = 1;
        end
        # define localizing matrices
        LocConst = Array{Any}(undef, 0); NumLoc = 0
        for i = 1:n
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "==";
            LocConst[NumLoc]["ord"] = 0;
        end
    elseif options["level"] == 0 && isequal(options["clique"], "off") # dense Shor's relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n));
        MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        for i = 1:n
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[i] = Dict();
            LocConst[i]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
            LocConst[i]["basis"] = sparse(I(n));
            LocConst[i]["typ"] = "==";
            LocConst[i]["ord"] = 0;
        end
    elseif options["level"] > 0 && isequal(options["clique"], "on") # sparse sublevel relaxation
        MomConst = Array{Any}(undef, 0); NumMom = 0;
        for i = 1:length(sizes)
            if options["level"] < sizes[i]
                NumMom += 1
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[cliques[i]', :]);
                MomConst[NumMom]["ord"] = options["ord"] - 1;
            else
                NumMom += 1
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[cliques[i]', :]);
                MomConst[NumMom]["ord"] = options["ord"];
            end
        end
        LocConst = Array{Any}(undef, 0); NumLoc = 0;
        for i = 1:n
            for j = 1:length(sizes)
                if i in cliques[j]
                    # x_i in clique j, find the index of x_i in clique j
                    idx_i = findall(cliques[j]'.==i)
                    clique = [cliques[j] cliques[j]];
                    if options["level"] < sizes[j]
                        # indices except x_i in clique j
                        set = filter(x->x∉idx_i, collect(1:sizes[j]))
                        set = [set; set]
                        if options["depth"] == 0
                            NumLoc += 1;
                            LocConst = vcat(LocConst, Array{Any}(undef, 1));
                            LocConst[NumLoc] = Dict();
                            LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                            LocConst[NumLoc]["basis"] = sparse(I(n));
                            LocConst[NumLoc]["ord"] = 0;
                            LocConst[NumLoc]["typ"] = "==";
                        else
                            depth = minimum([options["depth"] sizes[j]-1]);
                            subM = []; subL = [];
                            idx = heuristic(subM, subL, set, idx_i, depth, options; clique = cliques[j], cliques = cliques, siz = sizes[j], sizes = sizes)
                            for k = 1:depth
                                subset = vcat(idx_i, set[idx[k]:idx[k]+options["level"]-2])
                                NumMom += 1;
                                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                                MomConst[NumMom] = Dict();
                                MomConst[NumMom]["basis"] = sparse(E[sort(clique[subset]), :]);
                                MomConst[NumMom]["ord"] = options["ord"];
                                NumLoc += 1;
                                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                                LocConst[NumLoc] = Dict();
                                LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                                LocConst[NumLoc]["basis"] = sparse(E[sort(clique[subset]), :]);
                                LocConst[NumLoc]["ord"] = options["ord"] - 1;
                                LocConst[NumLoc]["typ"] = "==";
                            end
                        end
                    else
                        NumLoc += 1;
                        LocConst = vcat(LocConst, Array{Any}(undef, 1));
                        LocConst[NumLoc] = Dict();
                        LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                        LocConst[NumLoc]["basis"] = sparse(E[cliques[j]', :]);
                        LocConst[NumLoc]["ord"] = options["ord"] - 1;
                        LocConst[NumLoc]["typ"] = "==";
                    end
                    break
                end
            end
        end
    elseif options["level"] > 0 && isequal(options["clique"], "off") # dense sublevel relaxation
        MomConst = Array{Any}(undef, 1);
        LocConst = Array{Any}(undef, 0);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n));
        MomConst[1]["ord"] = 1;
        NumMom = 1; NumLoc = 0;
        clique = [1:n; 1:n]';
        if options["level"] < n
            for i = 1:n
                set = filter(x->x∉i, collect(1:n))
                set = [set; set]
                if options["depth"] == 0 # Shor's relaxation
                    NumLoc += 1;
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                    LocConst[NumLoc]["basis"] = sparse(I(n));
                    LocConst[NumLoc]["ord"] = 0;
                    LocConst[NumLoc]["typ"] = "==";
                else
                    depth = options["depth"];
                    idx = heuristic(M, L, set, i, depth, options)
                    for k = 1:depth
                        subset = vcat(i, set[idx[k]:idx[k]+options["level"]-2])
                        NumMom += 1;
                        MomConst = vcat(MomConst, Array{Any}(undef, 1));
                        MomConst[NumMom] = Dict();
                        MomConst[NumMom]["basis"] = sparse(E[sort(clique[subset]), :]);
                        MomConst[NumMom]["ord"] = options["ord"];
                        NumLoc += 1;
                        LocConst = vcat(LocConst, Array{Any}(undef, 1));
                        LocConst[NumLoc] = Dict();
                        LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                        LocConst[NumLoc]["basis"] = sparse(E[sort(clique[subset]), :]);
                        LocConst[NumLoc]["ord"] = options["ord"] - 1;
                        LocConst[NumLoc]["typ"] = "==";
                    end
                end
            end
        else
            NumMom += 1;
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[1:n, :]);
            MomConst[NumMom]["ord"] = options["ord"];
            for i = 1:n
                NumLoc += 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[1:n, :]);
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
            end
        end
    end
    @printf("\nMAX-CUT Problem: %d vertices, %d edges, clique %s, order %d, level %d, depth %d\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), uppercase(options["clique"]), options["ord"], options["level"], options["level"])
    OptVal, running_time, stat, M = solve_moment(typ, obj, MomConst, LocConst, options; cliques = cliques, sizes = sizes);
    return OptVal, running_time, stat, M
end

function solve_moment_mac(A, options)
    n = size(A, 1); typ = "max";
    E = sparse([I(n); I(n)]);
    obj = [zeros(1,n+1); [zeros(n,1) A]];
    if options["level"] == 0 # Shor's relaxation
        MomConst = Array{Any}(undef, 0); NumMom = 0
        NumMom += 1
        MomConst = vcat(MomConst, Array{Any}(undef, 1));
        MomConst[NumMom] = Dict();
        MomConst[NumMom]["basis"] = sparse(I(n));
        MomConst[NumMom]["ord"] = 1;
        LocConst = Array{Any}(undef, 0); NumLoc = 0
        for i = 1:n
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = 0;
        end
        NumLoc += 1
        LocConst = vcat(LocConst, Array{Any}(undef, 1));
        LocConst[NumLoc] = Dict();
        LocConst[NumLoc]["pol"] = sparse([-1 zeros(1,n); ones(n,1) E[1:n,:]]);
        LocConst[NumLoc]["basis"] = sparse(I(n));
        LocConst[NumLoc]["typ"] = "==";
        LocConst[NumLoc]["ord"] = 0;
    elseif options["level"] > 0
        MomConst = Array{Any}(undef, 1); LocConst = Array{Any}(undef, 0);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n));
        MomConst[1]["ord"] = 1;
        NumMom = 1; NumLoc = 0;
        clique = [1:n; 1:n]';
        if options["level"] < n && options["depth"] == 0
            for i = 1:n
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[i] = Dict();
                LocConst[i]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[i]["basis"] = sparse(I(n));
                LocConst[i]["typ"] = ">=";
                LocConst[i]["ord"] = 0;
            end
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[n+1] = Dict();
            LocConst[n+1]["pol"] = sparse([-1 zeros(1,n); ones(n,1) E[1:n,:]]);
            LocConst[n+1]["basis"] = sparse(I(n));
            LocConst[n+1]["typ"] = "==";
            LocConst[n+1]["ord"] = 0;
        elseif options["level"] < n && options["depth"] > 0
            depth = options["depth"];
            for i = 1:n
                set = filter(x->x∉i, collect(1:n))
                set = [set; set];
                idx = i;
                subset = vcat(i, set[idx[1]:idx[1]+options["level"]-2])
                NumMom = NumMom + 1;
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[sort(clique[subset]), :]);
                MomConst[NumMom]["ord"] = options["ord"];
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[sort(clique[subset]), :]);
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
                LocConst[NumLoc]["typ"] = ">=";
            end
            for k = 1:depth
                set = [1:n; 1:n]';
                idx = collect(1:1+depth-1);
                subset = set[idx[k]:idx[k]+options["level"]-1]
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([-1 zeros(1,n); ones(n,1) E[1:n,:]]);
                LocConst[NumLoc]["basis"] = sparse(E[sort(clique[subset]), :]);
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
                LocConst[NumLoc]["typ"] = "==";
            end
        else
            NumMom = NumMom + 1;
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[1:n, :]);
            MomConst[NumMom]["ord"] = options["ord"];
            for i = 1:n
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[1:n, :]);
                LocConst[NumLoc]["typ"] = ">=";
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
            end
            NumLoc = NumLoc + 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([-1 zeros(1,n); ones(n,1) E[1:n,:]]);
            LocConst[NumLoc]["basis"] = sparse(E[1:n, :]);
            LocConst[NumLoc]["typ"] = "==";
            LocConst[NumLoc]["ord"] = options["ord"] - 1;
        end
    end
    @printf("\nMAXIMUM-CLIQUE Problem: %d vertices, %d edges, order %d, level %d, depth %d\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), options["ord"], options["level"], options["depth"])
    OptVal, running_time, stat = solve_moment(typ, obj, MomConst, LocConst, options);
    return OptVal, running_time, stat
end

function solve_moment_mip(A, options; cliques = [], sizes = [], M = [])
    n = size(A, 1); typ = "min";
    E = sparse([I(n); I(n)]);
    obj = [zeros(1,n+1); [zeros(n,1) A]];
    if options["level"] == 0 && isequal(options["clique"], "on") # sparse Shor's relaxation
        # define moment matrices
        MomConst = Array{Any}(undef, 0); NumMom = 0
        for i = 1:length(sizes)
            NumMom += 1
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[cliques[i]', :]);
            MomConst[NumMom]["ord"] = 1;
        end
        # define localizing matrices
        LocConst = Array{Any}(undef, 0); NumLoc = 0
        for i = 1:n
            NumLoc += 1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "==";
            LocConst[NumLoc]["ord"] = 0;
        end
    elseif options["level"] == 0 && isequal(options["clique"], "off") # dense Shor's relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n));
        MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        for i = 1:n
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[i] = Dict();
            LocConst[i]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
            LocConst[i]["basis"] = sparse(I(n));
            LocConst[i]["typ"] = "==";
            LocConst[i]["ord"] = 0;
        end
    elseif options["level"] > 0 && isequal(options["clique"], "on") # sparse sublevel relaxation
        MomConst = Array{Any}(undef, 0); NumMom = 0;
        for i = 1:length(sizes)
            if options["level"] < sizes[i]
                NumMom += 1
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[cliques[i]', :]);
                MomConst[NumMom]["ord"] = options["ord"] - 1;
            else
                NumMom += 1
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[cliques[i]', :]);
                MomConst[NumMom]["ord"] = options["ord"];
            end
        end
        LocConst = Array{Any}(undef, 0); NumLoc = 0;
        for i = 1:n
            for j = 1:length(sizes)
                if i in cliques[j]
                    # x_i in clique j, find the index of x_i in clique j
                    idx_i = findall(cliques[j]'.==i)
                    clique = [cliques[j] cliques[j]];
                    if options["level"] < sizes[j]
                        # indices except x_i in clique j
                        set = filter(x->x∉idx_i, collect(1:sizes[j]))
                        set = [set; set]
                        if options["depth"] == 0
                            NumLoc += 1;
                            LocConst = vcat(LocConst, Array{Any}(undef, 1));
                            LocConst[NumLoc] = Dict();
                            LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                            LocConst[NumLoc]["basis"] = sparse(I(n));
                            LocConst[NumLoc]["ord"] = 0;
                            LocConst[NumLoc]["typ"] = "==";
                        else
                            depth = minimum([options["depth"] sizes[j]-1]);
                            subM = []; subL = [];
                            idx = heuristic(subM, subL, set, idx_i, depth, options; clique = cliques[j], cliques = cliques, siz = sizes[j], sizes = sizes)
                            for k = 1:depth
                                subset = vcat(idx_i, set[idx[k]:idx[k]+options["level"]-2])
                                NumMom += 1;
                                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                                MomConst[NumMom] = Dict();
                                MomConst[NumMom]["basis"] = sparse(E[sort(clique[subset]), :]);
                                MomConst[NumMom]["ord"] = options["ord"];
                                NumLoc += 1;
                                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                                LocConst[NumLoc] = Dict();
                                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                                LocConst[NumLoc]["basis"] = sparse(E[sort(clique[subset]), :]);
                                LocConst[NumLoc]["ord"] = options["ord"] - 1;
                                LocConst[NumLoc]["typ"] = "==";
                            end
                        end
                    else
                        NumLoc += 1;
                        LocConst = vcat(LocConst, Array{Any}(undef, 1));
                        LocConst[NumLoc] = Dict();
                        LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                        LocConst[NumLoc]["basis"] = sparse(E[cliques[j]', :]);
                        LocConst[NumLoc]["ord"] = options["ord"] - 1;
                        LocConst[NumLoc]["typ"] = "==";
                    end
                    break
                end
            end
        end
    elseif options["level"] > 0 && isequal(options["clique"], "off") # dense sublevel relaxation
        MomConst = Array{Any}(undef, 1);
        LocConst = Array{Any}(undef, 0);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n));
        MomConst[1]["ord"] = 1;
        NumMom = 1; NumLoc = 0;
        clique = [1:n; 1:n]';
        lv = options["level"];
        if options["level"] < n
            for i = 1:n
                set = filter(x->x∉i, collect(1:n))
                set = [set; set]
                if options["depth"] == 0 # Shor's relaxation
                    NumLoc += 1;
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                    LocConst[NumLoc]["basis"] = sparse(I(n));
                    LocConst[NumLoc]["ord"] = 0;
                    LocConst[NumLoc]["typ"] = "==";
                else
                    depth = options["depth"];
                    idx = heuristic(M, L, set, i, depth, options)
                    for k = 1:depth
                        subset = vcat(i, set[idx[k]:idx[k]+options["level"]-2])
                        NumMom += 1;
                        MomConst = vcat(MomConst, Array{Any}(undef, 1));
                        MomConst[NumMom] = Dict();
                        MomConst[NumMom]["basis"] = sparse(E[sort(clique[subset]), :]);
                        MomConst[NumMom]["ord"] = options["ord"];
                        NumLoc += 1;
                        LocConst = vcat(LocConst, Array{Any}(undef, 1));
                        LocConst[NumLoc] = Dict();
                        LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                        LocConst[NumLoc]["basis"] = sparse(E[sort(clique[subset]), :]);
                        LocConst[NumLoc]["ord"] = options["ord"] - 1;
                        LocConst[NumLoc]["typ"] = "==";
                    end
                end
            end
        else
            NumMom += 1;
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[1:n, :]);
            MomConst[NumMom]["ord"] = options["ord"];
            for i = 1:n
                NumLoc += 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[1:n, :]);
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
            end
        end
    end
    @printf("\nMIP Problem: %d vertices, %d edges, clique %s, order %d, level %d, depth %d\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), uppercase(options["clique"]), options["ord"], options["level"], options["depth"])
    OptVal, running_time, stat, M = solve_moment(typ, obj, MomConst, LocConst, options; cliques = cliques, sizes = sizes);
    return OptVal, running_time, stat, M
end

function solve_moment_qcqp(A, b, options)
    n = size(A, 1); typ = "min";
    E = sparse([I(n); I(n)]);
    obj = [0 b'; [b A]];
    if options["level"] == 0 && isequal(options["clique"], "off") # Shor's relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n));
        MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        NumMom = 1; NumLoc = 0;
        for i = 1:n
            NumLoc = NumLoc + 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = 0;
        end
        NumLoc = NumLoc + 1;
        LocConst = vcat(LocConst, Array{Any}(undef, 1));
        LocConst[NumLoc] = Dict();
        LocConst[NumLoc]["pol"] = sparse([-1 zeros(1,n); ones(n,1) 2*E[1:n,:]]);
        LocConst[NumLoc]["basis"] = sparse(I(n));
        LocConst[NumLoc]["typ"] = "==";
        LocConst[NumLoc]["ord"] = 0;
    elseif options["level"] > 0 && isequal(options["clique"], "off")
        MomConst = Array{Any}(undef, 1); LocConst = Array{Any}(undef, 0);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n));
        MomConst[1]["ord"] = 1;
        NumMom = 1; NumLoc = 0;
        clique = [1:n; 1:n]';
        if options["level"] < n && options["depth"] == 0
            for i = 1:n
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = ">=";
                LocConst[NumLoc]["ord"] = 0;
            end
            NumLoc = NumLoc + 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([-1 zeros(1,n); ones(n,1) 2*E[1:n,:]]);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "==";
            LocConst[NumLoc]["ord"] = 0;
        elseif options["level"] < n && options["depth"] > 0
            depth = options["depth"];
            for i = 1:n
                set = filter(x->x∉i, collect(1:n))
                set = [set; set];
                idx = i;
                subset = vcat(i, set[idx[1]:idx[1]+options["level"]-2])
                NumMom = NumMom + 1;
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[sort(clique[subset]), :]);
                MomConst[NumMom]["ord"] = options["ord"];
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[sort(clique[subset]), :]);
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
                LocConst[NumLoc]["typ"] = ">=";
            end
            for k = 1:depth
                set = [1:n; 1:n]';
                idx = collect(1:1+depth-1);
                subset = set[idx[k]:idx[k]+options["level"]-1]
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([-1 zeros(1,n); ones(n,1) 2*E[1:n,:]]);
                LocConst[NumLoc]["basis"] = sparse(E[sort(clique[subset]), :]);
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
                LocConst[NumLoc]["typ"] = "==";
            end
        else
            NumMom = NumMom + 1;
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(I(n));
            MomConst[NumMom]["ord"] = options["ord"];
            for i = 1:n
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = ">=";
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
            end
            NumLoc = NumLoc + 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([-1 zeros(1,n); ones(n,1) 2*E[1:n,:]]);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "==";
            LocConst[NumLoc]["ord"] = options["ord"] - 1;
        end
    end
    @printf("\nQCQP problem: %d variables, %d entries, clique %s, order %d, level %d, depth %d\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), uppercase(options["clique"]), options["ord"], options["level"], options["depth"])
    OptVal, running_time, stat = solve_moment(typ, obj, MomConst, LocConst, options);
    return OptVal, running_time, stat
end

# binary assignment, id: 0032
function solve_moment_qplib_BinaryAssignment(obj, constr, type, var,  options)
    typ = "min";
    n = length(var["u"]);
    E = sparse([I(n); I(n)]);
    objective = [obj["c"] obj["b"]'/2; [obj["b"]/2 obj["A"]]];
    if options["level"] == 0 && isequal(options["clique"], "off") # dense  Shor's relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n));
        MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        NumMom = 1; NumLoc = 0;
        for i = 1:length(constr["A"])
            idx = findall(!iszero, constr["b"][i][:]);
            pol = [constr["b"][i][idx] E[idx, :]];
            for j = 1:n
                for k = j:n
                    if constr["A"][i][j,k] != 0 && j == k
                        pol = vcat(pol, [constr["A"][i][j,k] 2*E[j,:]'])
                    elseif constr["A"][i][j,k] != 0 && j != k
                        pol = vcat(pol, [2*constr["A"][i][j,k] (E[j,:]+E[k,:])'])
                    end
                end
            end
            if constr["c"][i,1] != constr["c"][i,2] && constr["c"][i,1] == -Inf && constr["c"][i,2] != Inf # ∑y_i ≤ 5
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,2] zeros(1,n)]);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = "<=";
                LocConst[NumLoc]["ord"] = 0;
            elseif constr["c"][i,1] != constr["c"][i,2] && constr["c"][i,1] != -Inf && constr["c"][i,2] != Inf # -1 ≤ x_i - y_i ≤ 0
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([pol; 1 E[idx[1],:]'; -1 E[idx[2], :]'; 1 2*E[idx[1],:]'; 1 2*E[idx[2],:]'; -2 (E[idx[1],:].+E[idx[2],:])']);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = "<=";
                LocConst[NumLoc]["ord"] = 0;
            elseif constr["c"][i,1] == constr["c"][i,2] # ∑x_i = 1
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,1] zeros(1,n)]);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = 0;
            else
                error("NotImplementedError")
            end
        end
        for i = 1:n
            if var["type"][i] == 0 && var["l"][i] != -Inf && var["u"][i] != Inf # 0 ≤ x_i ≤ 1
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([var["l"][i]*var["u"][i] zeros(1,n); -var["l"][i]-var["u"][i] E[i,:]'; 1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = "<=";
                LocConst[NumLoc]["ord"] = 0;
            elseif var["type"][i] == 1 # y_i ∈ {0, 1}
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = 0;
            else
                error("NotImplementedError")
            end
        end
    elseif options["level"] > 0 && isequal(options["clique"], "off") # dense  sublevel relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n));
        MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        NumMom = 1; NumLoc = 0;
        for i = 1:length(constr["A"])
            idx = findall(!iszero, constr["b"][i][:]);
            pol = [constr["b"][i][idx] E[idx, :]];
            for j = 1:n
                for k = j:n
                    if constr["A"][i][j,k] != 0 && j == k
                        pol = vcat(pol, [constr["A"][i][j,k] 2*E[j,:]'])
                    elseif constr["A"][i][j,k] != 0 && j != k
                        pol = vcat(pol, [2*constr["A"][i][j,k] (E[j,:]+E[k,:])'])
                    end
                end
            end
            if constr["c"][i,1] != constr["c"][i,2] && constr["c"][i,1] == -Inf && constr["c"][i,2] != Inf # ∑y_i ≤ 5
                for k = 1:options["depth"]
                    subset = collect(50+k:50+k+options["level"]-1)
                    NumMom = NumMom + 1;
                    MomConst = vcat(MomConst, Array{Any}(undef, 1));
                    MomConst[NumMom] = Dict();
                    MomConst[NumMom]["basis"] = sparse(E[sort(subset), :]);
                    MomConst[NumMom]["ord"] = options["ord"];
                    NumLoc = NumLoc + 1;
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] =  sparse([pol; -constr["c"][i,2] zeros(1,n)]);
                    LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
                    LocConst[NumLoc]["typ"] = "<=";
                    LocConst[NumLoc]["ord"] = options["ord"]-1;
                end
            elseif constr["c"][i,1] != constr["c"][i,2] && constr["c"][i,1] != -Inf && constr["c"][i,2] != Inf # -1 ≤ x_i - y_i ≤ 0
                subset = vcat(collect(idx[1]:idx[1]+Int(options["level"]/2)-1), collect(idx[2]-Int(options["level"]/2)+1:idx[2]))
                NumMom = NumMom + 1;
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[sort(subset), :]);
                MomConst[NumMom]["ord"] = options["ord"];
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([pol; 1 E[idx[1],:]'; -1 E[idx[2], :]'; 1 2*E[idx[1],:]'; 1 2*E[idx[2],:]'; -2 (E[idx[1],:].+E[idx[2],:])']);
                LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
                LocConst[NumLoc]["typ"] = "<=";
                LocConst[NumLoc]["ord"] = options["ord"]-1;
            elseif constr["c"][i,1] == constr["c"][i,2] # ∑x_i = 1
                for k = 1:options["depth"]
                    subset = collect(k:k+options["level"]-1)
                    NumMom = NumMom + 1;
                    MomConst = vcat(MomConst, Array{Any}(undef, 1));
                    MomConst[NumMom] = Dict();
                    MomConst[NumMom]["basis"] = sparse(E[sort(subset), :]);
                    MomConst[NumMom]["ord"] = options["ord"];
                    NumLoc = NumLoc + 1;
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,1] zeros(1,n)]);
                    LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
                    LocConst[NumLoc]["typ"] = "==";
                    LocConst[NumLoc]["ord"] = options["ord"]-1;
                end
            else
                error("NotImplementedError")
            end
        end
        for i = 1:50
            subset = vcat(collect(i:i+Int(options["level"]/2)-1), collect(50+i-Int(options["level"]/2)+1:50+i))
            NumLoc = NumLoc + 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([var["l"][i]*var["u"][i] zeros(1,n); -var["l"][i]-var["u"][i] E[i,:]'; 1 2*E[i,:]']);
            LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
            LocConst[NumLoc]["typ"] = "<=";
            LocConst[NumLoc]["ord"] = options["ord"]-1;
            NumLoc = NumLoc + 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[50+i,:]'; -1 2*E[50+i,:]']);
            LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
            LocConst[NumLoc]["typ"] = "==";
            LocConst[NumLoc]["ord"] = options["ord"]-1;
        end
    else
        error("NotImplementedError")
    end
    @printf("\nQCQP problem: %d variables, %d entries, clique %s, order %d, level %d, depth %d\n", n, (sum(Matrix(1*(obj["A"].!=0))) - sum(diag(Matrix(1*(obj["A"].!=0)))))/2 + sum(diag(Matrix(1*(obj["A"].!=0)))), uppercase(options["clique"]), options["ord"], options["level"], options["depth"])
    OptVal, running_time, stat = solve_moment(typ, objective, MomConst, LocConst, options);
    return OptVal, running_time, stat
end

# qcqp, id: 1661, 1675, 1773, 1535, 1703,
function solve_moment_qplib_QCQP(obj, constr, type, var,  options)
    if type == "minimize"
        typ = "min";
    elseif type == "maximize"
        typ = "max";
    end
    n = length(var["u"]);
    E = sparse([I(n); I(n)]);
    objective = [obj["c"] obj["b"]'/2; [obj["b"]/2 obj["A"]]];
    if options["level"] == 0 && isequal(options["clique"], "off") # dense  Shor's relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n));
        MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        NumMom = 1; NumLoc = 0;
        for i = 1:length(constr["A"])
            idx = findall(!iszero, constr["b"][i][:]);
            pol = [constr["b"][i][idx] E[idx, :]];
            for j = 1:n
                for k = j:n
                    if constr["A"][i][j,k] != 0 && j == k
                        pol = vcat(pol, [constr["A"][i][j,k] 2*E[j,:]'])
                    elseif constr["A"][i][j,k] != 0 && j != k
                        pol = vcat(pol, [2*constr["A"][i][j,k] (E[j,:]+E[k,:])'])
                    end
                end
            end
            if constr["c"][i,1] != constr["c"][i,2] && constr["c"][i,1] == -Inf && constr["c"][i,2] != Inf
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,2] zeros(1,n)]);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = "<=";
                LocConst[NumLoc]["ord"] = 0;
            elseif constr["c"][i,1] != constr["c"][i,2] && constr["c"][i,1] != -Inf && constr["c"][i,2] == Inf
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,1] zeros(1,n)]);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = ">=";
                LocConst[NumLoc]["ord"] = 0;
            elseif constr["c"][i,1] != constr["c"][i,2] && constr["c"][i,1] != -Inf && constr["c"][i,2] != Inf
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,1] zeros(1,n)]);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = ">=";
                LocConst[NumLoc]["ord"] = 0;
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,2] zeros(1,n)]);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = "<=";
                LocConst[NumLoc]["ord"] = 0;
            elseif constr["c"][i,1] == constr["c"][i,2]
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,1] zeros(1,n)]);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = 0;
            else
                error("NotImplementedError")
            end
        end
        for i = 1:n
            # 0 ≤ x_i ≤ 1
            NumLoc = NumLoc + 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = 0;
        end
    elseif options["level"] > 0 && isequal(options["clique"], "off") # dense  sublevel relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n));
        MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        NumMom = 1; NumLoc = 0;
        for i = 1:length(constr["A"])
            idx = findall(!iszero, constr["b"][i][:]);
            pol = [constr["b"][i][idx] E[idx, :]];
            for j = 1:n
                for k = j:n
                    if constr["A"][i][j,k] != 0 && j == k
                        pol = vcat(pol, [constr["A"][i][j,k] 2*E[j,:]'])
                    elseif constr["A"][i][j,k] != 0 && j != k
                        pol = vcat(pol, [2*constr["A"][i][j,k] (E[j,:]+E[k,:])'])
                    end
                end
            end
            if sum(constr["A"][i]) == 0
                ord = 1; depth = options["depth"];
            else
                ord = 0; depth = 1;
            end
            for k = 1:options["depth"]
                subset = collect(k:k+options["level"]-1);
                if constr["c"][i,1] != constr["c"][i,2] && constr["c"][i,1] == -Inf && constr["c"][i,2] != Inf
                    NumLoc = NumLoc + 1;
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,2] zeros(1,n)]);
                    LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
                    LocConst[NumLoc]["typ"] = "<=";
                    LocConst[NumLoc]["ord"] = ord;
                elseif constr["c"][i,1] != constr["c"][i,2] && constr["c"][i,1] != -Inf && constr["c"][i,2] != Inf
                    NumLoc = NumLoc + 1;
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,1] zeros(1,n)]);
                    LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
                    LocConst[NumLoc]["typ"] = ">=";
                    LocConst[NumLoc]["ord"] = ord;
                    NumLoc = NumLoc + 1;
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,2] zeros(1,n)]);
                    LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
                    LocConst[NumLoc]["typ"] = "<=";
                    LocConst[NumLoc]["ord"] = ord;
                elseif constr["c"][i,1] == constr["c"][i,2]
                    NumLoc = NumLoc + 1;
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,1] zeros(1,n)]);
                    LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
                    LocConst[NumLoc]["typ"] = "==";
                    LocConst[NumLoc]["ord"] = ord;
                else
                    error("NotImplementedError")
                end
            end
        end
        for i = 1:n
            subset = collect(i:i+options["level"]-1);
            NumMom = NumMom + 1;
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[sort(subset), :]);
            MomConst[NumMom]["ord"] = options["ord"];
            # 0 ≤ x_i ≤ 1
            NumLoc = NumLoc + 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
            LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
            LocConst[NumLoc]["typ"] = ">=";
            LocConst[NumLoc]["ord"] = options["ord"]-1;
        end
    end
    @printf("\nQCQP problem: %d variables, %d entries, clique %s, order %d, level %d, depth %d\n", n, (sum(Matrix(1*(obj["A"].!=0))) - sum(diag(Matrix(1*(obj["A"].!=0)))))/2 + sum(diag(Matrix(1*(obj["A"].!=0)))), uppercase(options["clique"]), options["ord"], options["level"], options["depth"])
    OptVal, running_time, stat = solve_moment(typ, objective, MomConst, LocConst, options);
    return OptVal, running_time, stat
end

# binary, id: 0067, 0633, 2512, 3762, 5935, 5944, 10072, 10073, 10074
function solve_moment_qplib_binary(obj, constr, type, var,  options)
    if type == "minimize"
        typ = "min";
    elseif type == "maximize"
        typ = "max";
    end
    n = length(var["type"]);
    E = sparse([I(n); I(n)]);
    objective = [obj["c"] obj["b"]'/2; [obj["b"]/2 obj["A"]]];
    if options["level"] == 0 && isequal(options["clique"], "off") # dense  Shor's relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n));
        MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        NumMom = 1; NumLoc = 0;
        for i = 1:length(constr["A"])
            idx = findall(!iszero, constr["b"][i][:]);
            pol = [constr["b"][i][idx] E[idx, :]];
            for j = 1:n
                for k = j:n
                    if constr["A"][i][j,k] != 0 && j == k
                        pol = vcat(pol, [constr["A"][i][j,k] 2*E[j,:]'])
                    elseif constr["A"][i][j,k] != 0 && j != k
                        pol = vcat(pol, [2*constr["A"][i][j,k] (E[j,:]+E[k,:])'])
                    end
                end
            end
            if constr["c"][i,1] != constr["c"][i,2] && constr["c"][i,1] == -Inf && constr["c"][i,2] != Inf
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,2] zeros(1,n)]);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = "<=";
                LocConst[NumLoc]["ord"] = 0;
            elseif constr["c"][i,1] != constr["c"][i,2] && constr["c"][i,1] != -Inf && constr["c"][i,2] != Inf
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,1] zeros(1,n)]);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = ">=";
                LocConst[NumLoc]["ord"] = 0;
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,2] zeros(1,n)]);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = "<=";
                LocConst[NumLoc]["ord"] = 0;
            elseif constr["c"][i,1] == constr["c"][i,2]
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,1] zeros(1,n)]);
                LocConst[NumLoc]["basis"] = sparse(I(n));
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = 0;
                println(pol)
            else
                error("NotImplementedError")
            end
        end
        for i = 1:n
            # x_i ∈ {0, 1}
            NumLoc = NumLoc + 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
            LocConst[NumLoc]["basis"] = sparse(I(n));
            LocConst[NumLoc]["typ"] = "==";
            LocConst[NumLoc]["ord"] = 0;
        end
    elseif options["level"] > 0 && isequal(options["clique"], "off") # dense  sublevel relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n));
        MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        NumMom = 1; NumLoc = 0;
        for i = 1:length(constr["A"])
            idx = findall(!iszero, constr["b"][i][:]);
            pol = [constr["b"][i][idx] E[idx, :]];
            for j = 1:n
                for k = j:n
                    if constr["A"][i][j,k] != 0 && j == k
                        pol = vcat(pol, [constr["A"][i][j,k] 2*E[j,:]'])
                    elseif constr["A"][i][j,k] != 0 && j != k
                        pol = vcat(pol, [2*constr["A"][i][j,k] (E[j,:]+E[k,:])'])
                    end
                end
            end
            for k = 1:options["depth"]
                if length(idx) == n
                    subset = collect(k:k+options["level"]-1);
                elseif length(idx) >= options["level"]
                    subset = idx[1:options["level"]]
                elseif length(idx) < options["level"]
                    set = filter(x->x∉idx, collect(1:n))
                    subset = vcat(set[k:k+options["level"]-length(idx)-1], idx)
                end
                if constr["c"][i,1] != constr["c"][i,2] && constr["c"][i,1] == -Inf && constr["c"][i,2] != Inf
                    NumLoc = NumLoc + 1;
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,2] zeros(1,n)]);
                    LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
                    LocConst[NumLoc]["typ"] = "<=";
                    LocConst[NumLoc]["ord"] = options["ord"]-1;
                elseif constr["c"][i,1] != constr["c"][i,2] && constr["c"][i,1] != -Inf && constr["c"][i,2] != Inf
                    NumLoc = NumLoc + 1;
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,1] zeros(1,n)]);
                    LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
                    LocConst[NumLoc]["typ"] = ">=";
                    LocConst[NumLoc]["ord"] = options["ord"]-1;
                    NumLoc = NumLoc + 1;
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,2] zeros(1,n)]);
                    LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
                    LocConst[NumLoc]["typ"] = "<=";
                    LocConst[NumLoc]["ord"] = options["ord"]-1;
                elseif constr["c"][i,1] == constr["c"][i,2]
                    NumLoc = NumLoc + 1;
                    LocConst = vcat(LocConst, Array{Any}(undef, 1));
                    LocConst[NumLoc] = Dict();
                    LocConst[NumLoc]["pol"] = sparse([pol; -constr["c"][i,1] zeros(1,n)]);
                    LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
                    LocConst[NumLoc]["typ"] = "==";
                    LocConst[NumLoc]["ord"] = options["ord"]-1;
                else
                    error("NotImplementedError")
                end
            end
        end
        for i = 1:n
            subset = collect(i:i+options["level"]-1);
            NumMom = NumMom + 1;
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[NumMom] = Dict();
            MomConst[NumMom]["basis"] = sparse(E[sort(subset), :]);
            MomConst[NumMom]["ord"] = options["ord"];
            # x_i ∈ {0, 1}
            NumLoc = NumLoc + 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[NumLoc] = Dict();
            LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
            LocConst[NumLoc]["basis"] = sparse(E[sort(subset), :]);
            LocConst[NumLoc]["typ"] = "==";
            LocConst[NumLoc]["ord"] = options["ord"]-1;
        end
    end
    @printf("\nQCQP problem: %d variables, %d entries, clique %s, order %d, level %d, depth %d\n", n, (sum(Matrix(1*(obj["A"].!=0))) - sum(diag(Matrix(1*(obj["A"].!=0)))))/2 + sum(diag(Matrix(1*(obj["A"].!=0)))), uppercase(options["clique"]), options["ord"], options["level"], options["depth"])
    OptVal, running_time, stat = solve_moment(typ, objective, MomConst, LocConst, options);
    return OptVal, running_time, stat
end
