function solve_moment_lip_one_layer(A, b, c, x00, eps, options)
    if isequal(options["range"], "global")
        ϵ = 10;
    else
        ϵ = eps;
    end
    d0 = size(A[1], 2); d = size(A[1], 1); n = 2*d0+d; typ = "max";
    lx = x00 .- ϵ*ones(d0,1); ux = x00 .+ ϵ*ones(d0,1);
    lz = (A[1] .* (A[1] .> 0)) * lx .+ (A[1] .* (A[1] .< 0)) * ux .+ b[1];
    uz = (A[1] .* (A[1] .> 0)) * ux .+ (A[1] .* (A[1] .< 0)) * lx .+ b[1];
    ly = 1*(lz .> 0); uy = 1*(uz .> 0);
    E0 = sparse([zeros(1+n, n-1) [2 zeros(1,n); ones(n,1) I(n)]]);
    E = sparse([I(n); I(n)]);
    obj = sparse([zeros(n+1,1) E0[1:1+n, n+1:2*n]]);
    for i = 2:d0+1
        obj = vcat(obj, [zeros(n-i+2,1) E0[1:n-i+2, n-i+2:2*n-i+1]]);
    end
    for i = 2+d0:d0+d+1
        obj = vcat(obj, [[zeros(d0+d+2-i,1); c[i-d0-1]*A[1][i-d0-1,:]] E0[1:n-i+2, n-i+2:2*n-i+1]]);
    end
    for i = 2+d0+d:2*d0+d+1
        obj = vcat(obj, [zeros(n-i+2,1) E0[1:n-i+2, n-i+2:2*n-i+1]]);
    end
    if options["level"] == 0 # Shor's relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        for i = 1:d0
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[(i-1)*2+1] = Dict();
            LocConst[(i-1)*2+1]["pol"] = sparse([1 zeros(1,n); -1 2*E[d0+d+i,:]']);
            LocConst[(i-1)*2+1]["basis"] = sparse(I(n));
            LocConst[(i-1)*2+1]["typ"] = ">=";
            LocConst[(i-1)*2+1]["ord"] = 0;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[(i-1)*2+2] = Dict();
            LocConst[(i-1)*2+2]["pol"] = sparse([lx[i]*ux[i] zeros(1,n); -lx[i]-ux[i] E[i,:]'; 1 2*E[i,:]']);
            LocConst[(i-1)*2+2]["basis"] = sparse(I(n));
            LocConst[(i-1)*2+2]["typ"] = "<=";
            LocConst[(i-1)*2+2]["ord"] = 0;
        end
        for i = 1:d
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[2*d0+(i-1)*3+1] = Dict();
            LocConst[2*d0+(i-1)*3+1]["pol"] = sparse([1 E[d0+i,:]'; -1 2*E[d0+i,:]']);
            LocConst[2*d0+(i-1)*3+1]["basis"] = sparse(I(n));
            LocConst[2*d0+(i-1)*3+1]["typ"] = "==";
            LocConst[2*d0+(i-1)*3+1]["ord"] = 0;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[2*d0+(i-1)*3+2] = Dict();
            LocConst[2*d0+(i-1)*3+2]["pol"] = sparse([ly[i]*uy[i] zeros(1,n); -ly[i]-uy[i] E[d0+i,:]'; 1 2*E[d0+i,:]']);
            LocConst[2*d0+(i-1)*3+2]["basis"] = sparse(I(n));
            LocConst[2*d0+(i-1)*3+2]["typ"] = "<=";
            LocConst[2*d0+(i-1)*3+2]["ord"] = 0;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[2*d0+(i-1)*3+3] = Dict();
            LocConst[2*d0+(i-1)*3+3]["pol"] = sparse([-1/2*b[1][i] zeros(1,n); b[1][i] E[d0+i,:]'; [-1/2*A[1][i,:] E[1:d0, :]]; [A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1: 2*d0+d]]]);
            LocConst[2*d0+(i-1)*3+3]["basis"] = sparse(I(n));
            LocConst[2*d0+(i-1)*3+3]["typ"] = ">=";
            LocConst[2*d0+(i-1)*3+3]["ord"] = 0;
        end
    elseif options["level"] > 0
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        for i = 1:d0
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[1+i] = Dict();
            MomConst[1+i]["basis"] = sparse(E[[d0+1:d0+options["level"]; d0+d+i], :]);
            MomConst[1+i]["ord"] = options["ord"];
            # t^2<=1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[(i-1)*2+1] = Dict();
            LocConst[(i-1)*2+1]["pol"] = sparse([1 zeros(1,n); -1 2*E[d0+d+i,:]']);
            LocConst[(i-1)*2+1]["basis"] = sparse(E[[d0+1:d0+options["level"]; d0+d+i], :]);
            LocConst[(i-1)*2+1]["typ"] = ">=";
            LocConst[(i-1)*2+1]["ord"] = options["ord"] - 1;
            # (x-lx)(x-ux)<=0
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[(i-1)*2+2] = Dict();
            LocConst[(i-1)*2+2]["pol"] = sparse([lx[i]*ux[i] zeros(1,n); -lx[i]-ux[i] E[i,:]'; 1 2*E[i,:]']);
            LocConst[(i-1)*2+2]["basis"] = sparse(E[[i; d0+1:d0+options["level"]], :]);
            LocConst[(i-1)*2+2]["typ"] = "<=";
            LocConst[(i-1)*2+2]["ord"] = options["ord"] - 1;
        end
        for i = 1:d
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[1+d0+i] = Dict();
            MomConst[1+d0+i]["basis"] = sparse(E[[1:options["level"]; d0+i], :]);
            MomConst[1+d0+i]["ord"] = options["ord"];
            # y*(y-1)==0
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[2*d0+(i-1)*3+1] = Dict();
            LocConst[2*d0+(i-1)*3+1]["pol"] = sparse([1 E[d0+i,:]'; -1 2*E[d0+i,:]']);
            LocConst[2*d0+(i-1)*3+1]["basis"] = sparse(E[[1:options["level"]; d0+i], :]);
            LocConst[2*d0+(i-1)*3+1]["typ"] = "==";
            LocConst[2*d0+(i-1)*3+1]["ord"] = options["ord"] - 1;
            # (y-ly)(y-uy)<=0
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[2*d0+(i-1)*3+2] = Dict();
            LocConst[2*d0+(i-1)*3+2]["pol"] = sparse([ly[i]*uy[i] zeros(1,n); -ly[i]-uy[i] E[d0+i,:]'; 1 2*E[d0+i,:]']);
            LocConst[2*d0+(i-1)*3+2]["basis"] = sparse(E[[1:options["level"]; d0+i], :]);
            LocConst[2*d0+(i-1)*3+2]["typ"] = "<=";
            LocConst[2*d0+(i-1)*3+2]["ord"] = options["ord"] - 1;
            # (y-1/2)*(Ax+b)>=0
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[2*d0+(i-1)*3+3] = Dict();
            LocConst[2*d0+(i-1)*3+3]["pol"] = sparse([-1/2*b[1][i] zeros(1,n); b[1][i] E[d0+i,:]'; [-1/2*A[1][i,:] E[1:d0, :]]; [A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1: 2*d0+d]]]);
            LocConst[2*d0+(i-1)*3+3]["basis"] = sparse(E[[1:options["level"]; d0+i], :]);
            LocConst[2*d0+(i-1)*3+3]["typ"] = ">=";
            LocConst[2*d0+(i-1)*3+3]["ord"] = options["ord"] - 1;
        end
    end
    @printf("\n%s Lipschitz constant estimation problem: clique %s, order %d, level %d\n", uppercasefirst(options["range"]), uppercase(options["clique"]), options["ord"], options["level"])
    OptVal, running_time, stat = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    return OptVal, running_time, stat
end
# vars = matread("net_1_7.mat");
# A = vars["A"]; b = vars["b"]; c = vars["c"]; x00 = vars["x00"]; eps = 0.1;
# options = Dict();
# options["range"] = "global"; options["level"] = 5; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
# OptVal, running_time, status = solve_moment_lip_one_layer(A, b, c, x00, eps, options);

# n = 3;
# G = zeros(3*n,3*n)
# for i = 1:n
#     for j = n+1:2*n
#         G[i,j] = 1; G[j,i] = 1;
#     end
# end
# for i = n+1:2*n
#     for j = [1:n; 2*n+1:3*n]
#         G[i,j] = 1; G[j,i] = 1;
#     end
# end
# for i = 2*n+1:3*n
#     for j = n+1:2*n
#         G[i,j] = 1; G[j,i] = 1;
#     end
# end
# c,s,nG,g,ng = ChordalExtension(G, turn="on")

function solve_moment_lip_one_layer_tssos(A, b, c, x00, eps, options)
    if isequal(options["range"], "global")
        ϵ = 10;
    else
        ϵ = eps;
    end
    d0 = size(A[1],2); d = size(A[1],1);
    lx = x00.-ϵ*ones(d0,1); ux = x00.+ϵ*ones(d0,1);
    lz = (A[1] .* (A[1] .> 0)) * lx .+ (A[1] .* (A[1] .< 0)) * ux .+ b[1];
    uz = (A[1] .* (A[1] .> 0)) * ux .+ (A[1] .* (A[1] .< 0)) * lx .+ b[1];
    ly = 1*(lz .> 0); uy = 1*(uz .> 0);
    @polyvar t[1:d0] x[1:d0] y[1:d]
    pop = Array{Any}(undef, 1+2*d0+3*d)
    pop[1] = -(t'*A[1]'*diagm(y[:,1])*c)[1]
    for i = 1:d0
        # @printf("%d\n", i)
        pop[1+i] = 1-t[i]^2;
        pop[1+d0+i] = -(x[i]-lx[i])*(x[i]-ux[i]);
    end
    for i = 1:d
        # @printf("%d\n", i)
        pop[1+2*d0+i] = (y[i]-1/2)*(A[1][i,:]'*x+b[1][i]);
        pop[1+2*d0+d+i] = -(y[i]-ly[i])*(y[i]-uy[i]);
        pop[1+2*d0+2*d+i] = y[i]*(y[i]-1);
    end
    @printf("\n%s Lipschitz constant estimation problem: method %s\n", uppercasefirst(options["range"]), uppercase(options["method"]))
    running_time = Dict(); OptVal = Dict();
    start = time();
    if options["first"] == true
        OptValFirst, sol, data = blockcpop_first(pop, [t;x;y], 2, numeq=d, method=options["method"], QUIET=options["silent"]);
        running_time["first"] = time()-start; OptVal["first"] = -OptValFirst;
        if options["higher"] == true
            OptValHigher = 1; iter = 0;
            while OptValHigher > 0
                start = time();
                OptValHigher, sol, data = blockcpop_higher!(data, method=options["method"], QUIET=options["silent"]);
                running_time["higher$(iter)"] = time()-start;
                iter = iter + 1; OptVal["higher$(iter)"] = -OptValHigher;
            end
        end
    end
    if options["mix"] == true
        n = length(t) + length(x) + length(y);
        m = length(pop) - 1;
        numeq = d;
        coe=Array{Vector{Float64}}(undef, m+1);
        supp=Array{SparseMatrixCSC}(undef, m+1);
        for k = 1:m+1
            mon = monomials(pop[k]);
            coe[k] = coefficients(pop[k]);
            lt = length(mon);
            ssupp = zeros(UInt8,n,lt);
            for i = 1:lt
                for j = 1:d0
                    ssupp[j,i] = MultivariatePolynomials.degree(mon[i],t[j]);
                    ssupp[d0+j,i] = MultivariatePolynomials.degree(mon[i],x[j]);
                end
                for j = 1:d
                    ssupp[2*d0+j,i] = MultivariatePolynomials.degree(mon[i],y[j]);
                end
            end
            supp[k] = sparse(ssupp);
        end
        rlorder = ones(Int, m+1);
        rlorder[1] = 2;
        cliques, cql, cliquesize = clique_cdecomp(n, m, supp, rlorder, alg="amd");
        ts = 5;
        mclique,I,ncc,lmc,blocks,cl,blocksize,ub,sizes,ssupp,lt,fbasis,gbasis=get_cblocks_mix(rlorder,m,supp,cliques,cql,cliquesize,ts=ts,method=options["method"],chor_alg="amd");
        start = time();
        OptValMixFirst,supp1=blockcpop_mix(n,m,rlorder,supp,coe,cliques,cql,cliquesize,mclique,I,ncc,lmc,blocks,cl,blocksize,numeq=numeq,ts=ts,QUIET=options["silent"],solve=true);
        running_time["mix_first"] = time()-start;
        OptVal["mix_first"] = -OptValMixFirst;
        if options["mix_higher"] == true
            OptValMixHigher = 1; iter = 0;
            while OptValMixHigher > 0
                blocks,cl,blocksize,ub,sizes,status=get_chblocks_mix!(I,supp1,ssupp,lt,fbasis,gbasis,mclique,lmc,cliques,cliquesize,blocks,cl,blocksize,ub,sizes,method=options["method"],chor_alg="amd")
                start = time()
                OptValMixHigher,supp1=blockcpop_mix(n,m,rlorder,supp,coe,cliques,cql,cliquesize,mclique,I,ncc,lmc,blocks,cl,blocksize,numeq=numeq,ts=ts,QUIET=options["silent"],solve=true)
                running_time["mix_higher"] = time()-start;
                iter = iter + 1; OptVal["mix_higher$(iter)"] = -OptValMixHigher;
            end
        end
    end
    return OptVal, running_time
end
# vars = matread("net_1_7.mat");
# A = vars["A"]; b = vars["b"]; c = vars["c"]; x00 = vars["x00"]; eps = 0.1;
# options = Dict();
# options["range"] = "local"; options["method"] = "clique"; options["silent"] = false; options["first"] = true; options["higher"] = false; options["mix"] = false; options["mix_higher"] = true;
# OptVal, running_time = solve_moment_lip_one_layer_tssos(A, b, c, x00, eps, options);

function solve_moment_lip_two_layer(A, b, c, x00, eps, options)
end
# options = Dict();
# options["range"] = "global"; options["level"] = 3; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = false;

function solve_moment_lip_two_layer_tssos(A, b, c, x00, eps, options)
    if isequal(options["range"], "global")
        ϵ = 10;
    else
        ϵ = eps;
    end
    d0 = size(A[1],2); d1 = size(A[1],1); d2 = size(A[2],1);
    lx = x00.-ϵ*ones(d0,1); ux = x00.+ϵ*ones(d0,1);
    lz1 = (A[1] .* (A[1] .> 0)) * lx .+ (A[1] .* (A[1] .< 0)) * ux .+ b[1];
    uz1 = (A[1] .* (A[1] .> 0)) * ux .+ (A[1] .* (A[1] .< 0)) * lx .+ b[1];
    ly1 = 1*(lz1 .> 0); uy1 = 1*(uz1 .> 0);
    lx1 = lz1.*(lz1.>=0); ux1 = uz1.*(uz1.>=0);
    lz2 = (A[2] .* (A[2] .> 0)) * lx1 .+ (A[2] .* (A[2] .< 0)) * ux1 .+ b[2];
    uz2 = (A[2] .* (A[2] .> 0)) * ux1 .+ (A[2] .* (A[2] .< 0)) * lx1 .+ b[2];
    ly2 = 1*(lz2 .> 0); uy2 = 1*(uz2 .> 0);
    @polyvar t[1:d0] x[1:d0] x1[1:d1] y1[1:d1] y2[1:d2]
    pop = Array{Any}(undef, 1+2*d0+7*d1+3*d2)
    pop[1] = -(t'*A[1]'*diagm(y1[:,1])*A[2]'*diagm(y2[:,1])*c)[1]
    for i = 1:d0
        # @printf("%d\n", i)
        pop[1+i] = 1-t[i]^2;
        pop[1+d0+i] = -(x[i]-lx[i])*(x[i]-ux[i]);
    end
    for i = 1:d1
        # @printf("%d\n", i)
        pop[1+2*d0+i] = (y1[i]-1/2)*(A[1][i,:]'*x+b[1][i]);
        pop[1+2*d0+d1+i] = -(y1[i]-ly1[i])*(y1[i]-uy1[i]);
        pop[1+2*d0+2*d1+i] = x1[i]-A[1][i,:]'*x-b[1][i];
        pop[1+2*d0+3*d1+i] = x1[i];
        pop[1+2*d0+4*d1+i] = -(x1[i]-lx1[i])*(x1[i]-ux1[i]);
    end
    for i = 1:d2
        pop[1+2*d0+5*d1+i] = (y2[i]-1/2)*(A[2][i,:]'*x1+b[2][i]);
        pop[1+2*d0+5*d1+d2+i] = -(y2[i]-ly2[i])*(y2[i]-uy2[i]);
        pop[1+2*d0+5*d1+2*d2+i] = y2[i]*(y2[i]-1);
    end
    for i = 1:d1
        pop[1+2*d0+5*d1+3*d2+i] = x1[i]*(x1[i]-A[1][i,:]'*x-b[1][i]);
        pop[1+2*d0+6*d1+3*d2+i] = y1[i]*(y1[i]-1);
    end
    @printf("\n%s Lipschitz constant estimation problem: method %s\n", uppercasefirst(options["range"]), uppercase(options["method"]))
    running_time = Dict(); OptVal = Dict();
    start = time();
    if options["first"] == true
        OptValFirst, sol, data = blockcpop_first(pop, [t;x;x1;y1;y2], 2, numeq=2*d1+d2, method=options["method"], QUIET=options["silent"]);
        running_time["first"] = time()-start; OptVal["first"] = -OptValFirst;
        if options["higher"] == true
            OptValHigher = 1; iter = 0;
            while OptValHigher > 0
                start = time();
                OptValHigher, sol, data = blockcpop_higher!(data, method=options["method"], QUIET=options["silent"]);
                running_time["higher$(iter)"] = time()-start;
                iter = iter + 1; OptVal["higher$(iter)"] = -OptValHigher;
            end
        end
    end
    if options["mix"] == true
        n = length(t) + length(x) + length(x1) + length(y1) + length(y2);
        m = length(pop) - 1;
        numeq = 2*d1+d2;
        coe=Array{Vector{Float64}}(undef, m+1);
        supp=Array{SparseMatrixCSC}(undef, m+1);
        for k = 1:m+1
            mon = monomials(pop[k]);
            coe[k] = coefficients(pop[k]);
            lt = length(mon);
            ssupp = zeros(UInt8,n,lt);
            for i = 1:lt
                for j = 1:d0
                    ssupp[j,i] = MultivariatePolynomials.degree(mon[i],t[j]);
                    ssupp[d0+j,i] = MultivariatePolynomials.degree(mon[i],x[j]);
                end
                for j = 1:d1
                    ssupp[2*d0+j,i] = MultivariatePolynomials.degree(mon[i],x1[j]);
                    ssupp[2*d0+d1+j,i] = MultivariatePolynomials.degree(mon[i],y1[j]);
                end
                for j = 1:d2
                    ssupp[2*d0+2*d1+j,i] = MultivariatePolynomials.degree(mon[i],y2[j]);
                end
            end
            supp[k] = sparse(ssupp);
        end
        rlorder = ones(Int, m+1);
        rlorder[1] = 2;
        cliques, cql, cliquesize = clique_cdecomp(n, m, supp, rlorder, alg="amd");
        ts = 5;
        mclique,I,ncc,lmc,blocks,cl,blocksize,ub,sizes,ssupp,lt,fbasis,gbasis=get_cblocks_mix(rlorder,m,supp,cliques,cql,cliquesize,ts=ts,method=options["method"],chor_alg="amd");
        start = time();
        OptValMixFirst,supp1=blockcpop_mix(n,m,rlorder,supp,coe,cliques,cql,cliquesize,mclique,I,ncc,lmc,blocks,cl,blocksize,numeq=numeq,ts=ts,QUIET=options["silent"],solve=true);
        running_time["mix_first"] = time()-start;
        OptVal["mix_first"] = -OptValMixFirst;
        if options["mix_higher"] == true
            OptValMixHigher = 1; iter = 0;
            while OptValMixHigher > 0
                blocks,cl,blocksize,ub,sizes,status=get_chblocks_mix!(I,supp1,ssupp,lt,fbasis,gbasis,mclique,lmc,cliques,cliquesize,blocks,cl,blocksize,ub,sizes,method=options["method"],chor_alg="amd")
                start = time()
                OptValMixHigher,supp1=blockcpop_mix(n,m,rlorder,supp,coe,cliques,cql,cliquesize,mclique,I,ncc,lmc,blocks,cl,blocksize,numeq=numeq,ts=ts,QUIET=options["silent"],solve=true)
                running_time["mix_higher"] = time()-start;
                iter = iter + 1; OptVal["mix_higher$(iter)"] = -OptValMixHigher;
            end
        end
    end
    return OptVal, running_time
end
# vars = matread("net_2_20.mat");
# A = vars["A"]; b = vars["b"]; c = vars["c"]; x00 = vars["x00"]; eps = 0.1;
# options = Dict();
# options["range"] = "local"; options["method"] = "clique"; options["silent"] = false; options["first"] = true; options["higher"] = false; options["mix"] = false; options["mix_higher"] = true;
# OptVal, running_time = solve_moment_lip_two_layer_tssos(A, b, c, x00, eps, options);

function solve_moment_cert(A, b, c, x00, eps, options)
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
    E0 = sparse([zeros(1+n, n-1) [2 zeros(1,n); ones(n,1) I(n)]]);
    E = sparse([I(n); I(n)]);
    obj = sparse([[zeros(d0+1,1); c] E0[1:1+n, n+1:2*n]]);
    for i = 2:n+1
        obj = vcat(obj, [zeros(n-i+2,1) E0[1:n-i+2, n-i+2:2*n-i+1]]);
    end
    if options["level"] == 0 # Shor's relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        for i = 1:d0
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[i] = Dict();
            LocConst[i]["pol"] = sparse([lx[i]*ux[i] zeros(1,n); -lx[i]-ux[i] E[i,:]'; 1 2*E[i,:]']);
            LocConst[i]["basis"] = sparse(I(n));
            LocConst[i]["typ"] = "<=";
            LocConst[i]["ord"] = 0;
        end
        for i = 1:d
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+1] = Dict();
            LocConst[d0+(i-1)*4+1]["pol"] = sparse([lx1[i]*ux1[i] zeros(1,n); -lx1[i]-ux1[i] E[d0+i,:]'; 1 2*E[d0+i,:]']);
            LocConst[d0+(i-1)*4+1]["basis"] = sparse(I(n));
            LocConst[d0+(i-1)*4+1]["typ"] = "<=";
            LocConst[d0+(i-1)*4+1]["ord"] = 0;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+2] = Dict();
            LocConst[d0+(i-1)*4+2]["pol"] = sparse([1 E[d0+i,:]']);
            LocConst[d0+(i-1)*4+2]["basis"] = sparse(I(n));
            LocConst[d0+(i-1)*4+2]["typ"] = ">=";
            LocConst[d0+(i-1)*4+2]["ord"] = 0;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+3] = Dict();
            LocConst[d0+(i-1)*4+3]["pol"] = sparse([1 E[d0+i,:]'; -A[1][i,:] E[1:d0, :]; -b[1][i] zeros(1,n)]);
            LocConst[d0+(i-1)*4+3]["basis"] = sparse(I(n));
            LocConst[d0+(i-1)*4+3]["typ"] = ">=";
            LocConst[d0+(i-1)*4+3]["ord"] = 0;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+4] = Dict();
            LocConst[d0+(i-1)*4+4]["pol"] = sparse([1 2*E[d0+i,:]'; [-A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1:d0+d]]; -b[1][i] E[d0+i,:]']);
            LocConst[d0+(i-1)*4+4]["basis"] = sparse(I(n));
            LocConst[d0+(i-1)*4+4]["typ"] = "==";
            LocConst[d0+(i-1)*4+4]["ord"] = 0;
        end
    elseif options["level"] > 0
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        for i = 1:d0
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[1+i] = Dict();
            MomConst[1+i]["basis"] = sparse(E[[i; d0+1:d0+options["level"]], :]);
            MomConst[1+i]["ord"] = options["ord"];
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[i] = Dict();
            LocConst[i]["pol"] = sparse([lx[i]*ux[i] zeros(1,n); -lx[i]-ux[i] E[i,:]'; 1 2*E[i,:]']);
            LocConst[i]["basis"] = sparse(E[[i; d0+1:d0+options["level"]], :]);
            LocConst[i]["typ"] = "<=";
            LocConst[i]["ord"] = options["ord"] - 1;
        end
        for i = 1:d
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[1+d0+i] = Dict();
            MomConst[1+d0+i]["basis"] = sparse(E[[1:options["level"]; d0+i], :]);
            MomConst[1+d0+i]["ord"] = options["ord"];
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+1] = Dict();
            LocConst[d0+(i-1)*4+1]["pol"] = sparse([lx1[i]*ux1[i] zeros(1,n); -lx1[i]-ux1[i] E[d0+i,:]'; 1 2*E[d0+i,:]']);
            LocConst[d0+(i-1)*4+1]["basis"] = sparse(E[[1:options["level"]; d0+i], :]);
            LocConst[d0+(i-1)*4+1]["typ"] = "<=";
            LocConst[d0+(i-1)*4+1]["ord"] = options["ord"] - 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+2] = Dict();
            LocConst[d0+(i-1)*4+2]["pol"] = sparse([1 E[d0+i,:]']);
            LocConst[d0+(i-1)*4+2]["basis"] = sparse(E[[1:options["level"]; d0+i], :]);
            LocConst[d0+(i-1)*4+2]["typ"] = ">=";
            LocConst[d0+(i-1)*4+2]["ord"] = options["ord"] - 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+3] = Dict();
            LocConst[d0+(i-1)*4+3]["pol"] = sparse([1 E[d0+i,:]'; -A[1][i,:] E[1:d0, :]; -b[1][i] zeros(1,n)]);
            LocConst[d0+(i-1)*4+3]["basis"] = sparse(E[[1:options["level"]; d0+i], :]);
            LocConst[d0+(i-1)*4+3]["typ"] = ">=";
            LocConst[d0+(i-1)*4+3]["ord"] = options["ord"] - 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+4] = Dict();
            LocConst[d0+(i-1)*4+4]["pol"] = sparse([1 2*E[d0+i,:]'; [-A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1:d0+d]]; -b[1][i] E[d0+i,:]']);
            LocConst[d0+(i-1)*4+4]["basis"] = sparse(E[[1:options["level"]; d0+i], :]);
            LocConst[d0+(i-1)*4+4]["typ"] = "==";
            LocConst[d0+(i-1)*4+4]["ord"] = options["ord"] - 1;
        end
    end
    @printf("\n%s Robustness certification problem: clique %s, order %d, level %d\n", uppercasefirst(options["range"]), uppercase(options["clique"]), options["ord"], options["level"])
    OptVal, running_time, stat = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    return OptVal, running_time, stat
end
# vars = matread("net_1_7.mat");
# # A = vars["A"]; b = vars["b"]; c = vars["c"]; x00 = vars["x00"]; eps = 0.1;
# A = Array{Any}(undef, 1); b = Array{Any}(undef, 1);
# A[1] = vars["A"]; b[1] = vars["b"]; c = vars["c"]; x00 = vars["x00"]; eps = 0.1;
# options = Dict();
# options["range"] = "local"; options["level"] = 0; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
# OptVal, running_time, status = solve_moment_cert(A, b, c, x00, eps, options);

# n = 3;
# G = zeros(2*n,2*n)
# for i = 1:n
#     for j = n+1:2*n
#         G[i,j] = 1; G[j,i] = 1;
#     end
# end
# c,s,nG,g,ng = ChordalExtension(G, turn="on")

function solve_moment_cert_tssos(A, b, c, x00, eps, options)
    if isequal(options["range"], "global")
        ϵ = 10;
    else
        ϵ = eps;
    end
    d0 = size(A[1],2); d = size(A[1],1);
    lx = x00.-ϵ*ones(d0,1); ux = x00.+ϵ*ones(d0,1);
    lz1 = (A[1] .* (A[1] .> 0)) * lx .+ (A[1] .* (A[1] .< 0)) * ux .+ b[1];
    uz1 = (A[1] .* (A[1] .> 0)) * ux .+ (A[1] .* (A[1] .< 0)) * lx .+ b[1];
    lx1 = lz1.*(lz1.>=0); ux1 = uz1.*(uz1.>=0);
    @polyvar x[1:d0] x1[1:d]
    pop = Array{Any}(undef, 1+d0+4*d)
    pop[1] = -(c'*x1)[1];
    for i = 1:d0
        # @printf("%d\n", i)
        pop[1+i] = -(x[i]-lx[i])*(x[i]-ux[i]);
    end
    for i = 1:d
        # @printf("%d\n", i)
        pop[1+d0+i] = -(x1[i]-lx1[i])*(x1[i]-ux1[i]);
        pop[1+d0+d+i] = x1[i];
        pop[1+d0+2*d+i] = x1[i]-A[1][i,:]'*x-b[1][i];
        pop[1+d0+3*d+i] = x1[i]*(x1[i]-A[1][i,:]'*x-b[1][i]);
    end
    @printf("\n%s Lipschitz constant estimation problem: method %s\n", uppercasefirst(options["range"]), uppercase(options["method"]))
    running_time = Dict(); OptVal = Dict();
    start = time();
    if options["first"] == true
        OptValFirst, sol, data = blockcpop_first(pop, [x;x1], 2, numeq=d, method=options["method"], QUIET=options["silent"]);
        running_time["first"] = time()-start; OptVal["first"] = -OptValFirst;
        if options["higher"] == true
            OptValHigher = 1; iter = 0;
            while OptValHigher > 0
                start = time();
                OptValHigher, sol, data = blockcpop_higher!(data, method=options["method"], QUIET=options["silent"]);
                running_time["higher$(iter)"] = time()-start;
                iter = iter + 1; OptVal["higher$(iter)"] = -OptValHigher;
            end
        end
    end
    if options["mix"] == true
        n = length(x) + length(x1);
        m = length(pop) - 1;
        numeq = d;
        coe=Array{Vector{Float64}}(undef, m+1);
        supp=Array{SparseMatrixCSC}(undef, m+1);
        for k = 1:m+1
            mon = monomials(pop[k]);
            coe[k] = coefficients(pop[k]);
            lt = length(mon);
            ssupp = zeros(UInt8,n,lt);
            for i = 1:lt
                for j = 1:d0
                    ssupp[j,i] = MultivariatePolynomials.degree(mon[i],x[j]);
                end
                for j = 1:d
                    ssupp[d0+j,i] = MultivariatePolynomials.degree(mon[i],x1[j]);
                end
            end
            supp[k] = sparse(ssupp);
        end
        rlorder = ones(Int, m+1);
        rlorder[1] = 2;
        cliques, cql, cliquesize = clique_cdecomp(n, m, supp, rlorder, alg="amd");
        ts = 5;
        mclique,I,ncc,lmc,blocks,cl,blocksize,ub,sizes,ssupp,lt,fbasis,gbasis=get_cblocks_mix(rlorder,m,supp,cliques,cql,cliquesize,ts=ts,method=options["method"],chor_alg="amd");
        start = time();
        OptValMixFirst,supp1=blockcpop_mix(n,m,rlorder,supp,coe,cliques,cql,cliquesize,mclique,I,ncc,lmc,blocks,cl,blocksize,numeq=numeq,ts=ts,QUIET=options["silent"],solve=true);
        running_time["mix_first"] = time()-start;
        OptVal["mix_first"] = -OptValMixFirst;
        if options["mix_higher"] == true
            OptValMixHigher = 1; iter = 0;
            while OptValMixHigher > 0
                blocks,cl,blocksize,ub,sizes,status=get_chblocks_mix!(I,supp1,ssupp,lt,fbasis,gbasis,mclique,lmc,cliques,cliquesize,blocks,cl,blocksize,ub,sizes,method=options["method"],chor_alg="amd")
                start = time()
                OptValMixHigher,supp1=blockcpop_mix(n,m,rlorder,supp,coe,cliques,cql,cliquesize,mclique,I,ncc,lmc,blocks,cl,blocksize,numeq=numeq,ts=ts,QUIET=options["silent"],solve=true)
                running_time["mix_higher"] = time()-start;
                iter = iter + 1; OptVal["mix_higher$(iter)"] = -OptValMixHigher;
            end
        end
    end
    return OptVal, running_time
end
# vars = matread("net_1_7.mat");
# A = vars["A"]; b = vars["b"]; c = vars["c"]; x00 = vars["x00"]; eps = 0.1;
# options = Dict();
# options["range"] = "local"; options["method"] = "clique"; options["silent"] = false; options["first"] = false; options["higher"] = false; options["mix"] = true; options["mix_higher"] = true;
# OptVal, running_time = solve_moment_cert_tssos(A, b, c, x00, eps, options);

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
        cliques, sizes = ChordalExtension(A);
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
                    MomConst[NumMom]["basis"] = sparse(E[clique[j:j+options["level"]-1], :]);
                    MomConst[NumMom]["ord"] = options["ord"];
                end
            else
                NumMom = NumMom + 1;
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[cliques[i], :]');
                MomConst[NumMom]["ord"] = options["ord"];
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
                    LocConst[NumLoc]["basis"] = sparse(E[clique[j:j+options["level"]-1], :]);
                    LocConst[NumLoc]["typ"] = "==";
                    LocConst[NumLoc]["ord"] = options["ord"] - 1;
                end
            else
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[cliques[idx], :]');
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
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
                MomConst[NumMom]["basis"] = sparse(E[clique[i:i+options["level"]-1], :]);
                MomConst[NumMom]["ord"] = options["ord"];
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[clique[i:i+options["level"]-1], :]);
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
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
                LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[1:n, :]);
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
            end
        end
    end
    @printf("\nMAX-CUT Problem: %d vertices, %d edges, clique %s, order %d, level %d\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), uppercase(options["clique"]), options["ord"], options["level"])
    OptVal, running_time, stat = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    return OptVal, running_time, stat
end
# vars = matread("maxcut.mat");
# A = vars["A"]; W = ones(size(A, 1), size(A, 1));
# options = Dict();
# options["level"] = 5; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
# OptVal, running_time, status = solve_moment_maxcut(A, W, options)

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
        cliques, sizes = ChordalExtension(A);
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
                    MomConst[NumMom]["basis"] = sparse(E[clique[j:j+options["level"]-1], :]);
                    MomConst[NumMom]["ord"] = options["ord"];
                end
            else
                NumMom = NumMom + 1;
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[cliques[i], :]');
                MomConst[NumMom]["ord"] = options["ord"];
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
                    LocConst[NumLoc]["basis"] = sparse(E[clique[j:j+options["level"]-1], :]);
                    LocConst[NumLoc]["typ"] = "==";
                    LocConst[NumLoc]["ord"] = options["ord"] - 1;
                end
            else
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[cliques[idx], :]');
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
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
                MomConst[NumMom]["basis"] = sparse(E[clique[i:i+options["level"]-1], :]);
                MomConst[NumMom]["ord"] = options["ord"];
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[clique[i:i+options["level"]-1], :]);
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
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
                LocConst[NumLoc]["typ"] = "==";
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
            end
        end
    end
    @printf("\nMixed Integer Programming (MIP): %d variables, %d entries, clique %s, order %d, level %d\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), uppercase(options["clique"]), options["ord"], options["level"])
    OptVal, running_time, stat = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    return OptVal, running_time, stat
end
# vars = matread("mip.mat");
# A = vars["L"];
# options = Dict();
# options["level"] = 15; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
# OptVal, running_time, status = solve_moment_mip(A, options);
