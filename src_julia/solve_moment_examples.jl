function solve_moment_lip_one_layer(A, b, c, x00, e, options)
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
    # println(lx); println(ux); println(lx.*ux)
    E0 = sparse([zeros(1+n, n-1) [2 zeros(1,n); ones(n,1) I(n)]]);
    E = sparse([I(n); I(n)]);
    obj = [zeros(1+d0,n+1); [zeros(d,1+d0+d) 1/2*diagm(c[:])*A[1]; zeros(d0,1+d0) 1/2*A[1]'*diagm(c[:]) zeros(d0,d0)]];
    # obj = sparse([zeros(n+1,1) E0[1:1+n, n+1:2*n]]);
    # for i = 2:d0+1
    #     obj = vcat(obj, [zeros(n-i+2,1) E0[1:n-i+2, n-i+2:2*n-i+1]]);
    # end
    # for i = 2+d0:d0+d+1
    #     obj = vcat(obj, [[zeros(d0+d+2-i,1); c[i-d0-1]*A[1][i-d0-1,:]] E0[1:n-i+2, n-i+2:2*n-i+1]]);
    # end
    # for i = 2+d0+d:2*d0+d+1
    #     obj = vcat(obj, [zeros(n-i+2,1) E0[1:n-i+2, n-i+2:2*n-i+1]]);
    # end
    if options["level"] == 0 # Shor's relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        for i = 1:d0
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[(i-1)*2+1] = Dict();
            LocConst[(i-1)*2+1]["pol"] = sparse([1 zeros(1,n); -1 2*E[d0+d+i,:]']); #t^2<=1
            LocConst[(i-1)*2+1]["basis"] = sparse(I(n));
            LocConst[(i-1)*2+1]["typ"] = ">=";
            LocConst[(i-1)*2+1]["ord"] = 0;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[(i-1)*2+2] = Dict();
            # LocConst[(i-1)*2+2]["pol"] = sparse([lx[i]*ux[i]*100 zeros(1,n); -(lx[i]+ux[i])*10 E[i,:]'; 1*1 2*E[i,:]']); #(x-lx)(x-ux)<=0
            # LocConst[(i-1)*2+2]["pol"] = sparse([-1 E[i,:]'; 1 2*E[i,:]']); #x(x-1)<=0
            LocConst[(i-1)*2+2]["pol"] = sparse([-1 zeros(1,n); 1 2*E[i,:]']); #x^2<=1
            LocConst[(i-1)*2+2]["basis"] = sparse(I(n));
            LocConst[(i-1)*2+2]["typ"] = "<=";
            LocConst[(i-1)*2+2]["ord"] = 0;
        end
        for i = 1:d
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[2*d0+(i-1)*3+1] = Dict();
            LocConst[2*d0+(i-1)*3+1]["pol"] = sparse([1 E[d0+i,:]'; -1 2*E[d0+i,:]']); #y*(y-1)==0
            LocConst[2*d0+(i-1)*3+1]["basis"] = sparse(I(n));
            LocConst[2*d0+(i-1)*3+1]["typ"] = "==";
            LocConst[2*d0+(i-1)*3+1]["ord"] = 0;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[2*d0+(i-1)*3+2] = Dict();
            LocConst[2*d0+(i-1)*3+2]["pol"] = sparse([ly[i]*uy[i] zeros(1,n); -(ly[i]+uy[i]) E[d0+i,:]'; 1 2*E[d0+i,:]']); #(y-ly)(y-uy)<=0
            LocConst[2*d0+(i-1)*3+2]["basis"] = sparse(I(n));
            LocConst[2*d0+(i-1)*3+2]["typ"] = "<=";
            LocConst[2*d0+(i-1)*3+2]["ord"] = 0;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[2*d0+(i-1)*3+3] = Dict();
            # LocConst[2*d0+(i-1)*3+3]["pol"] = sparse([-1/2*b[1][i]*10 zeros(1,n); b[1][i]*10 E[d0+i,:]'; [-1/2*A[1][i,:] E[1:d0, :]]; [A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1: 2*d0+d]]]); #(y-1/2)*(Ax+b)>=0
            # LocConst[2*d0+(i-1)*3+3]["pol"] = sparse([-1/2*(b[1][i]+(A[1][i,:]'*lx)[1]) zeros(1,n); b[1][i]+(A[1][i,:]'*lx)[1] E[d0+i,:]'; -1/2*A[1][i,:].*(ux-lx) E[1:d0, :]; A[1][i,:].*(ux-lx) E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1: 2*d0+d]]); #(y-1/2)*(Ax(ux-lx)+A*lx+b)>=0
            LocConst[2*d0+(i-1)*3+3]["pol"] = sparse([-1/2*(b[1][i]+(A[1][i,:]'*(ux+lx)/2)[1]) zeros(1,n); b[1][i]+(A[1][i,:]'*(ux+lx)/2)[1] E[d0+i,:]'; -1/2*A[1][i,:].*(ux-lx)/2 E[1:d0, :]; A[1][i,:].*(ux-lx)/2 E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1: 2*d0+d]]); #(y-1/2)*(Ax(ux-lx)/2+A*(ux+lx)/2+b)>=0
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
            if i+options["level"]-1<=n
                cl = clique[i:i+options["level"]-1]
            else
                cl = clique[[1:i+options["level"]-1-n; i:n]]
            end
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[1+(i-1)*2+1] = Dict();
            MomConst[1+(i-1)*2+1]["basis"] = sparse(E[unique([d0+1:d0+options["level"]; d0+d+i]), :]);
            # MomConst[1+(i-1)*2+1]["basis"] = sparse(E[d0+d+i:d0+d+i+1, :]);
            MomConst[1+(i-1)*2+1]["ord"] = options["ord"];
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[1+(i-1)*2+2] = Dict();
            MomConst[1+(i-1)*2+2]["basis"] = sparse(E[unique([i; d0+1:d0+options["level"]]), :]);
            # MomConst[1+(i-1)*2+2]["basis"] = sparse(E[i:i+1, :]);
            MomConst[1+(i-1)*2+2]["ord"] = options["ord"];
            # t^2<=1
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[(i-1)*2+1] = Dict();
            LocConst[(i-1)*2+1]["pol"] = sparse([1 zeros(1,n); -1 2*E[d0+d+i,:]']); #t^2<=1
            # println(sparse([1 zeros(1,n); -1 2*E[d0+d+i,:]']))
            LocConst[(i-1)*2+1]["basis"] = sparse(E[unique([d0+1:d0+options["level"]; d0+d+i]), :]);
            # LocConst[(i-1)*2+1]["basis"] = sparse(E[d0+d+i:d0+d+i+1, :]);
            LocConst[(i-1)*2+1]["typ"] = ">=";
            LocConst[(i-1)*2+1]["ord"] = options["ord"] - 1;
            # (x-lx)(x-ux)<=0
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[(i-1)*2+2] = Dict();
            # LocConst[(i-1)*2+2]["pol"] = sparse([lx[i]*ux[i] zeros(1,n); -(lx[i]+ux[i]) E[i,:]'; 1 2*E[i,:]']); #(x-lx)(x-ux)<=0
            # LocConst[(i-1)*2+2]["pol"] = sparse([-1 E[i,:]'; 1 2*E[i,:]']); #x(x-1)<=0
            LocConst[(i-1)*2+2]["pol"] = sparse([-1 zeros(1,n); 1 2*E[i,:]']); #x^2<=1
            # println(sparse([lx[i]*ux[i] zeros(1,n); -lx[i]-ux[i] E[i,:]'; 1 2*E[i,:]']))
            LocConst[(i-1)*2+2]["basis"] = sparse(E[unique([i; d0+1:d0+options["level"]]), :]);
            # LocConst[(i-1)*2+2]["basis"] = sparse(E[i:i+1, :]);
            LocConst[(i-1)*2+2]["typ"] = "<=";
            LocConst[(i-1)*2+2]["ord"] = options["ord"] - 1;
        end
        for i = 1:d
            MomConst = vcat(MomConst, Array{Any}(undef, 1));
            MomConst[1+2*d0+i] = Dict();
            MomConst[1+2*d0+i]["basis"] = sparse(E[unique([1:options["level"]; d0+i]), :]);
            # MomConst[1+2*d0+i]["basis"] = sparse(E[d0+i:d0+i+options["level"]-1, :]);
            MomConst[1+2*d0+i]["ord"] = options["ord"];
            # y*(y-1)==0
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[2*d0+(i-1)*3+1] = Dict();
            LocConst[2*d0+(i-1)*3+1]["pol"] = sparse([1 E[d0+i,:]'; -1 2*E[d0+i,:]']); #y*(y-1)==0
            # println(sparse([1 E[d0+i,:]'; -1 2*E[d0+i,:]']))
            LocConst[2*d0+(i-1)*3+1]["basis"] = sparse(E[unique([1:options["level"]; d0+i]), :]);
            # LocConst[2*d0+(i-1)*3+1]["basis"] = sparse(E[d0+i:d0+i+1, :]);
            LocConst[2*d0+(i-1)*3+1]["typ"] = "==";
            LocConst[2*d0+(i-1)*3+1]["ord"] = options["ord"] - 1;
            # (y-ly)(y-uy)<=0
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[2*d0+(i-1)*3+2] = Dict();
            LocConst[2*d0+(i-1)*3+2]["pol"] = sparse([ly[i]*uy[i] zeros(1,n); -(ly[i]+uy[i]) E[d0+i,:]'; 1 2*E[d0+i,:]']); #(y-ly)(y-uy)<=0
            # println(sparse([ly[i]*uy[i] zeros(1,n); -ly[i]-uy[i] E[d0+i,:]'; 1 2*E[d0+i,:]']))
            LocConst[2*d0+(i-1)*3+2]["basis"] = sparse(E[unique([1:options["level"]; d0+i]), :]);
            # LocConst[2*d0+(i-1)*3+2]["basis"] = sparse(E[d0+i:d0+i+1, :]);
            LocConst[2*d0+(i-1)*3+2]["typ"] = "<=";
            LocConst[2*d0+(i-1)*3+2]["ord"] = options["ord"] - 1;
            # (y-1/2)*(Ax+b)>=0
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[2*d0+(i-1)*3+3] = Dict();
            # LocConst[2*d0+(i-1)*3+3]["pol"] = sparse([-1/2*b[1][i] zeros(1,n); b[1][i] E[d0+i,:]'; [-1/2*A[1][i,:] E[1:d0, :]]; [A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1: 2*d0+d]]]); #(y-1/2)*(Ax+b)>=0
            # LocConst[2*d0+(i-1)*3+3]["pol"] = sparse([-1/2*(b[1][i]+(A[1][i,:]'*lx)[1]) zeros(1,n); b[1][i]+(A[1][i,:]'*lx)[1] E[d0+i,:]'; -1/2*A[1][i,:].*(ux-lx) E[1:d0, :]; A[1][i,:].*(ux-lx) E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1: 2*d0+d]]); #(y-1/2)*(Ax(ux-lx)+A*lx+b)>=0
            LocConst[2*d0+(i-1)*3+3]["pol"] = sparse([-1/2*(b[1][i]+(A[1][i,:]'*(ux+lx)/2)[1]) zeros(1,n); b[1][i]+(A[1][i,:]'*(ux+lx)/2)[1] E[d0+i,:]'; -1/2*A[1][i,:].*(ux-lx)/2 E[1:d0, :]; A[1][i,:].*(ux-lx)/2 E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1: 2*d0+d]]); #(y-1/2)*(Ax(ux-lx)/2+A*(ux+lx)/2+b)>=0
            # println(sparse([-1/2*b[1][i] zeros(1,n); b[1][i] E[d0+i,:]'; [-1/2*A[1][i,:] E[1:d0, :]]; [A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1: 2*d0+d]]]))
            LocConst[2*d0+(i-1)*3+3]["basis"] = sparse(E[unique([1:options["level"]; d0+i]), :]);
            # LocConst[2*d0+(i-1)*3+3]["basis"] = sparse(E[d0+i:d0+i+options["level"]-1, :]);
            LocConst[2*d0+(i-1)*3+3]["typ"] = ">=";
            LocConst[2*d0+(i-1)*3+3]["ord"] = options["ord"] - 1;
        end
    end
    @printf("\n%s Lipschitz constant estimation problem: clique %s, order %d, level %d\n", uppercasefirst(options["range"]), uppercase(options["clique"]), options["ord"], options["level"])
    # solve_moment_manual(typ, obj, MomConst, LocConst, options);
    OptVal, running_time, stat = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    return OptVal, running_time, stat
end
# vars = matread("net_1_15.mat"); A = vars["A"]; b = vars["b"]; c = vars["c"]; x00 = vars["x00"]; ϵ = 0.1; options = Dict(); options["range"] = "local"; options["level"] = 6; options["clique"] = "off"; options["ord"] = 2; options["silent"] = false; options["quad"] = true; OptVal, running_time, status = solve_moment_lip_one_layer(A, b, c, x00, ϵ, options);

# vars["MomMat"][2][1:3,1:3] == 1*vars["MomMat"][1][[1,4,6],[1,4,6]]
# vars["MomMat"][3][1:3,1:3] == 1*vars["MomMat"][1][[1,2,4],[1,2,4]]
# vars["MomMat"][4][1:3,1:3] == 1*vars["MomMat"][1][[1,4,7],[1,4,7]]
# vars["MomMat"][5][1:3,1:3] == 1*vars["MomMat"][1][[1,3,4],[1,3,4]]
# vars["MomMat"][6][1:3,1:3] == 1*vars["MomMat"][1][[1,2,4],[1,2,4]]
# vars["MomMat"][7][1:3,1:3] == 1*vars["MomMat"][1][[1,2,5],[1,2,5]]
# vars["MomMat"][3] == vars["MomMat"][6]
# vars["MomMat"][3][[1,3,6],[1,3,6]] == vars["MomMat"][2][[1,2,4],[1,2,4]]
# vars["MomMat"][3][[1,3,6],[1,3,6]] == vars["MomMat"][6][[1,3,6],[1,3,6]]
# vars["MomMat"][3][[1,3,6],[1,3,6]] == vars["MomMat"][5][[1,3,6],[1,3,6]]
# vars["MomMat"][2][[1,2,4],[1,2,4]] == vars["MomMat"][4][[1,2,4],[1,2,4]]
# vars["MomMat"][3][[1,2,4],[1,2,4]] == vars["MomMat"][6][[1,2,4],[1,2,4]]
# vars["MomMat"][3][[1,2,4],[1,2,4]] == vars["MomMat"][7][[1,2,4],[1,2,4]]

# benchmarks = ["net_1_5.mat", "net_1_10.mat", "net_1_15.mat", "net_1_20.mat", "net_1_25.mat", "net_1_30.mat", "net_1_35.mat", "net_1_40.mat"]#, "net_1_45.mat", "net_1_50.mat", "net_1_55.mat", "net_1_60.mat", "net_1_65.mat", "net_1_70.mat", "net_1_75.mat", "net_1_80.mat", "net_1_85.mat", "net_1_90.mat", "net_1_95.mat", "net_1_100.mat"]
# obj = Dict(); t = Dict(); options = Dict(); options["range"] = "local"; options["clique"] = "off"; options["ord"] = 2; options["silent"] = false; options["quad"] = true;
# for ben in benchmarks
#     vars = matread(ben);
#     A = Array{Any}(undef, 1); b = Array{Any}(undef, 1);
#     A[1] = vars["A"][1]; b[1] = vars["b"][1]; c = vars["c"]; x00 = vars["x00"]; ϵ = 0.1;
#     for l in [4 8]
#         options["level"] = l;
#         # OptVal, running_time, status = solve_moment_cert(A, b, c, x00, ϵ, options);
#         OptVal, running_time, status = solve_moment_lip_one_layer(A, b, c, x00, ϵ, options);
#         obj[@sprintf("%s_level_%d", ben, l)] = OptVal; t[@sprintf("%s_level_%d", ben, l)] = running_time["solv"];
#     end
# end

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
# vars = matread("net_1_5.mat");
# A = vars["A"]; b = vars["b"]; c = vars["c"]; x00 = vars["x00"]; ϵ = 0.1;
# options = Dict();
# options["range"] = "global"; options["method"] = "clique"; options["silent"] = false; options["first"] = false; options["higher"] = true; options["mix"] = true; options["mix_higher"] = false;
# OptVal, running_time = solve_moment_lip_one_layer_tssos(A, b, c, x00, ϵ, options);

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
# A = vars["A"]; b = vars["b"]; c = vars["c"]; x00 = vars["x00"]; ϵ = 0.1;
# options = Dict();
# options["range"] = "local"; options["method"] = "clique"; options["silent"] = false; options["first"] = true; options["higher"] = false; options["mix"] = false; options["mix_higher"] = true;
# OptVal, running_time = solve_moment_lip_two_layer_tssos(A, b, c, x00, ϵ, options);

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
    obj = [zeros(1,d0+1) 1/2*c'; zeros(d0,1+n); 1/2*c zeros(d,n)];
    # obj = sparse([[zeros(d0+1,1); c] E0[1:1+n, n+1:2*n]]);
    # for i = 2:n+1
    #     obj = vcat(obj, [zeros(n-i+2,1) E0[1:n-i+2, n-i+2:2*n-i+1]]);
    # end
    if options["level"] == 0 # Shor's relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
        LocConst = Array{Any}(undef, 0);
        for i = 1:d0
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[i] = Dict();
            LocConst[i]["pol"] = sparse([lx[i]*ux[i] zeros(1,n); -lx[i]-ux[i] E[i,:]'; 1 2*E[i,:]']); #(x-lx)(x-ux)<=0
            LocConst[i]["basis"] = sparse(I(n));
            LocConst[i]["typ"] = "<=";
            LocConst[i]["ord"] = 0;
        end
        for i = 1:d
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+1] = Dict();
            LocConst[d0+(i-1)*4+1]["pol"] = sparse([lx1[i]*ux1[i] zeros(1,n); -lx1[i]-ux1[i] E[d0+i,:]'; 1 2*E[d0+i,:]']); #(x1-lx1)(x1-ux1)<=0
            LocConst[d0+(i-1)*4+1]["basis"] = sparse(I(n));
            LocConst[d0+(i-1)*4+1]["typ"] = "<=";
            LocConst[d0+(i-1)*4+1]["ord"] = 0;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+2] = Dict();
            LocConst[d0+(i-1)*4+2]["pol"] = sparse([1 E[d0+i,:]']); #x1>=0
            LocConst[d0+(i-1)*4+2]["basis"] = sparse(I(n));
            LocConst[d0+(i-1)*4+2]["typ"] = ">=";
            LocConst[d0+(i-1)*4+2]["ord"] = 0;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+3] = Dict();
            LocConst[d0+(i-1)*4+3]["pol"] = sparse([1 E[d0+i,:]'; -A[1][i,:] E[1:d0, :]; -b[1][i] zeros(1,n)]); #x1>=Ax+b
            LocConst[d0+(i-1)*4+3]["basis"] = sparse(I(n));
            LocConst[d0+(i-1)*4+3]["typ"] = ">=";
            LocConst[d0+(i-1)*4+3]["ord"] = 0;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+4] = Dict();
            LocConst[d0+(i-1)*4+4]["pol"] = sparse([1 2*E[d0+i,:]'; [-A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1:d0+d]]; -b[1][i] E[d0+i,:]']); #x1(x1-Ax-b)=0
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
            LocConst[i]["pol"] = sparse([lx[i]*ux[i] zeros(1,n); -lx[i]-ux[i] E[i,:]'; 1 2*E[i,:]']); #(x-lx)(x-ux)<=0
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
            LocConst[d0+(i-1)*4+1]["pol"] = sparse([lx1[i]*ux1[i] zeros(1,n); -lx1[i]-ux1[i] E[d0+i,:]'; 1 2*E[d0+i,:]']); #(x1-lx1)(x1-ux1)<=0
            LocConst[d0+(i-1)*4+1]["basis"] = sparse(E[[1:options["level"]; d0+i], :]);
            LocConst[d0+(i-1)*4+1]["typ"] = "<=";
            LocConst[d0+(i-1)*4+1]["ord"] = options["ord"] - 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+2] = Dict();
            LocConst[d0+(i-1)*4+2]["pol"] = sparse([1 E[d0+i,:]']); #x1>=0
            LocConst[d0+(i-1)*4+2]["basis"] = sparse(E[[1:options["level"]; d0+i], :]);
            LocConst[d0+(i-1)*4+2]["typ"] = ">=";
            LocConst[d0+(i-1)*4+2]["ord"] = options["ord"] - 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+3] = Dict();
            LocConst[d0+(i-1)*4+3]["pol"] = sparse([1 E[d0+i,:]'; -A[1][i,:] E[1:d0, :]; -b[1][i] zeros(1,n)]); #x1>=Ax+b
            LocConst[d0+(i-1)*4+3]["basis"] = sparse(E[[1:options["level"]; d0+i], :]);
            LocConst[d0+(i-1)*4+3]["typ"] = ">=";
            LocConst[d0+(i-1)*4+3]["ord"] = options["ord"] - 1;
            LocConst = vcat(LocConst, Array{Any}(undef, 1));
            LocConst[d0+(i-1)*4+4] = Dict();
            LocConst[d0+(i-1)*4+4]["pol"] = sparse([1 2*E[d0+i,:]'; [-A[1][i,:] E[1:d0, 1:d0+i-1] ones(d0, 1) E[1:d0, d0+i+1:d0+d]]; -b[1][i] E[d0+i,:]']); #x1(x1-Ax-b)=0
            LocConst[d0+(i-1)*4+4]["basis"] = sparse(E[[1:options["level"]; d0+i], :]);
            LocConst[d0+(i-1)*4+4]["typ"] = "==";
            LocConst[d0+(i-1)*4+4]["ord"] = options["ord"] - 1;
        end
    end
    @printf("\n%s Robustness certification problem: clique %s, order %d, level %d\n", uppercasefirst(options["range"]), uppercase(options["clique"]), options["ord"], options["level"])
    OptVal, running_time, stat = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    return OptVal, running_time, stat
end
# vars = matread("net_1_5.mat"); A = Array{Any}(undef, 1); b = Array{Any}(undef, 1); A[1] = vars["A"][1]; b[1] = vars["b"][1]; c = vars["c"]; x00 = vars["x00"]; ϵ = 0.1; options = Dict(); options["range"] = "local"; options["level"] = 2; options["clique"] = "off"; options["ord"] = 2; options["silent"] = false; options["quad"] = true; OptVal, running_time, status = solve_moment_cert(A, b, c, x00, ϵ, options);

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
    @printf("\n%s Robustness certification problem: method %s\n", uppercasefirst(options["range"]), uppercase(options["method"]))
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
# vars = matread("net_1_5.mat");
# A = vars["A"]; b = vars["b"]; c = vars["c"]; x00 = vars["x00"]; ϵ = 0.1;
# options = Dict();
# options["range"] = "global"; options["method"] = "clique"; options["silent"] = false; options["first"] = false; options["higher"] = false; options["mix"] = true; options["mix_higher"] = true;
# OptVal, running_time = solve_moment_cert_tssos(A, b, c, x00, ϵ, options);

function solve_moment_maxcut(A, W, options)
    n = size(A, 1); typ = "max";
    L = sparse(diagm(((A.*W)*ones(n,1))[:,1])) - A.*W;
    E0 = sparse([zeros(1+n, n-1) [2 zeros(1,n); ones(n,1) I(n)]]);
    E = sparse([I(n); I(n)]);
    # L = L/2;
    # for i = 1:n
    #     L[i,i] = L[i,i]/2;
    # end
    # obj = sparse([zeros(n+1,1) E0[1:1+n, n+1:2*n]]);
    # for i = 2:n+1
    #     obj = vcat(obj, [L[i-1:n,i-1] E0[1:n-i+2, n-i+2:2*n-i+1]]);
    # end
    obj = [zeros(1,n+1); [zeros(n,1) 1/4*L]];
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
        # lv = 2; set = unique(rand(1:n,50))
        lv = options["level"]; set = 1:n
        if options["level"] < n
            for i = 1:n
                NumMom = NumMom + 1;
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                if i in set
                    MomConst[NumMom]["basis"] = sparse(E[clique[i:i+options["level"]-1], :]);
                    MomConst[NumMom]["ord"] = options["ord"];
                else
                    MomConst[NumMom]["basis"] = sparse(E[clique[i:i+lv-1], :]);
                    MomConst[NumMom]["ord"] = options["ord"];
                end
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                if i in set
                    LocConst[NumLoc]["basis"] = sparse(E[clique[i:i+options["level"]-1], :]);
                    LocConst[NumLoc]["ord"] = options["ord"] - 1;
                else
                    LocConst[NumLoc]["basis"] = sparse(E[clique[i:i+lv-1], :]);
                    LocConst[NumLoc]["ord"] = options["ord"] - 1;
                end
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
# vars = matread("maxcut_g20.txt.mat");
# A = vars["A"]; W = vars["W"];
# options = Dict();
# options["level"] = 6; options["clique"] = "off"; options["ord"] = 2; options["silent"] = false; options["quad"] = true; OptVal, running_time, status = solve_moment_maxcut(A, W, options)

function solve_moment_maxcut_multi_order(A, W, options)
    n = size(A, 1); typ = "max";
    L = sparse(diagm(((A.*W)*ones(n,1))[:,1])) - A.*W;
    E0 = sparse([zeros(1+n, n-1) [2 zeros(1,n); ones(n,1) I(n)]]);
    E = sparse([I(n); I(n)]);
    obj = [zeros(1,n+1); [zeros(n,1) 1/4*L]];
    iter = 0;
    # Shor's relaxation
    MomConst = Array{Any}(undef, 1);
    MomConst[1] = Dict();
    MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
    LocConst = Array{Any}(undef, 0);
    clique = [1:n; 1:n]';
    for i = 1:n
        MomConst = vcat(MomConst, Array{Any}(undef, 1));
        MomConst[i+1] = Dict();
        MomConst[i+1]["basis"] = sparse(E[1:n, :]);
        MomConst[i+1]["ord"] = 0;
        LocConst = vcat(LocConst, Array{Any}(undef, 1));
        LocConst[i] = Dict();
        LocConst[i]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
        LocConst[i]["basis"] = sparse(I(n));
        LocConst[i]["typ"] = "==";
        LocConst[i]["ord"] = 0;
    end
    @printf("\nMAX-CUT Problem: %d vertices, %d edges, level %d, order %d, epoch %d\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), options["level"], options["ord"], iter)
    OptVal, running_time, stat, M = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    deg = A*ones(n,1)[:,1];
    while iter <= options["iter"]
        # F = svd(M[2:end, 2:end]); U = F.U; u = U[:,1];
        # idx = findall(abs.(u.^2 .- 1) .== maximum(abs.(u.^2 .- 1)))
        idx = findall(deg .== maximum(deg));
        if isempty(idx)
            break
        else
            MaxNum = minimum([50, length(idx)])
            @printf("Targeting vertex: "); println(idx[1:MaxNum]); deleteat!(deg, idx[1:MaxNum])
            for i in idx[1:MaxNum]
                MomConst[i+1]["basis"] = sparse(E[clique[i:i+options["level"]-1], :]);
                MomConst[i+1]["ord"] = options["ord"];
                LocConst[i]["pol"] = sparse([1 zeros(1,n); -1 2*E[i,:]']);
                LocConst[i]["basis"] = sparse(E[clique[i:i+options["level"]-1], :]);
                LocConst[i]["typ"] = "==";
                LocConst[i]["ord"] = options["ord"] - 1;
            end
            @printf("\nMAX-CUT Problem: %d vertices, %d edges, level %d, order %d, epoch %d\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), options["level"], options["ord"], iter+1)
            OptVal, running_time, stat, M = solve_moment_manual(typ, obj, MomConst, LocConst, options);
        end
        iter = iter + 1;
    end
    @printf("\nOperation terminated at epoch %d, optimal value is %.2f, solver status is: ", iter-1, OptVal); println(stat)
    return OptVal, stat, iter
end
# vars = matread("maxcut_G11.mat");
# A = vars["A"]; W = vars["W"];
# options = Dict();
# options["iter"] = 10; options["level"] = 6; options["ord"] = 2; options["silent"] = true; options["quad"] = true; OptVal, running_time, status = solve_moment_maxcut_multi_order(A, W, options)

function solve_moment_maxcut_tssos(A, W, options)
    d = size(A,1)
    L = sparse(diagm(((A.*W)*ones(d,1))[:,1])) - A.*W
    @polyvar x[1:d]
    pop = Array{Any}(undef, 1+d)
    pop[1] = -1/4*x'*L*x;
    for i = 1:d
        # @printf("%d\n", i)
        pop[i+1] = x[i]^2-1
    end
    @printf("\nMAX-CUT Problem: %d vertices, %d edges, method %s\n", d, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), uppercase(options["method"]))
    running_time = Dict(); OptVal = Dict();
    start = time();
    if options["first"] == true
        OptValFirst, sol, data = blockcpop_first(pop, x, 2, numeq=d, method=options["method"], QUIET=options["silent"]);
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
        n = length(x);
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
                for j = 1:d
                    ssupp[j,i] = MultivariatePolynomials.degree(mon[i],x[j]);
                end
            end
            supp[k] = sparse(ssupp);
        end
        rlorder = ones(Int, m+1);
        rlorder[1] = 2;
        cliques, cql, cliquesize = clique_cdecomp(n, m, supp, rlorder, alg="amd");
        ts = 5;
        mclique,I,ncc,lmc,blocks,cl,blocksize,ub,sizes,ssupp,lt,fbasis,gbasis=get_cblocks_mix(rlorder,m,supp,cliques,cql,cliquesize,ts=ts,method=options["method"],chor_alg="max");
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
# vars = matread("maxcut_G11.mat");
# A = vars["A"]; W = vars["W"];
# options = Dict();
# options["method"] = "block"; options["silent"] = false; options["first"] = false; options["higher"] = true; options["mix"] = true; options["mix_higher"] = true; OptVal, running_time = solve_moment_maxcut_tssos(A, W, options)

# n=6
# @polyvar x[1:n]
# f=1+sum(x.^4)+x[1]*x[2]*x[3]+x[3]*x[4]*x[5]+x[3]*x[4]*x[6]+x[3]*x[5]*x[6]+x[4]*x[5]*x[6]
# pop=[f,1-sum(x[1:3].^2),1-sum(x[1:4].^2)]
# m=length(pop)-1
# coe=Array{Vector{Float64}}(undef, m+1)
# supp=Array{SparseMatrixCSC}(undef, m+1)
# for k=1:m+1
#     mon=monomials(pop[k])
#     coe[k]=coefficients(pop[k])
#     lt=length(mon)
#     ssupp=zeros(UInt8,n,lt)
#     for i=1:lt
#         for j=1:n
#             ssupp[j,i]=MultivariatePolynomials.degree(mon[i],x[j])
#         end
#     end
#     supp[k]=sparse(ssupp)
# end
# dg=2*ones(Int,m)
# order=2 # the relaxation order of Lasserre hierarchy
# opt,sol,data=cs_tssos_first(n,m,dg,supp,coe,order,numeq=0,TS="block")
# opt,sol,data=cs_tssos_higher!(data,TS="block")


#
# benchmarks = ["maxcut_g05_60.0.mat", "maxcut_g05_80.0.mat", "maxcut_g05_100.0.mat", "maxcut_pm1d_80.0.mat", "maxcut_pm1d_100.0.mat", "maxcut_pm1s_80.0.mat", "maxcut_pm1s_100.0.mat", "maxcut_pw01_100.0.mat", "maxcut_pw05_100.0.mat", "maxcut_pw09_100.0.mat", "maxcut_w01_100.0.mat", "maxcut_w05_100.0.mat", "maxcut_w09_100.0.mat"]
# for ben in benchmarks
#     for l = [0 4 6 8]
#         vars = matread(ben)
#         A = vars["A"]; W = vars["W"];
#         options = Dict();
#         options["level"] = l; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
#         OptVal, running_time, status = solve_moment_maxcut(A, W, options)
#     end
# end

function solve_moment_mip(A, options)
    n = size(A, 1); typ = "min";
    E0 = sparse([zeros(1+n, n-1) [2 zeros(1,n); ones(n,1) I(n)]]);
    E = sparse([I(n); I(n)]);
    # A = 2*A;
    # for i = 1:n
    #     A[i,i] = A[i,i]/2;
    # end
    # obj = sparse([zeros(n+1,1) E0[1:1+n, n+1:2*n]]);
    # for i = 2:n+1
    #     obj = vcat(obj, [A[i-1:n,i-1] E0[1:n-i+2, n-i+2:2*n-i+1]]);
    # end
    obj = [zeros(1,n+1); [zeros(n,1) A]];
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
# vars = matread("mip_gka10b.sparse.mat");
# A = vars["L"];
# options = Dict();
# options["level"] = 0; options["clique"] = "off"; options["ord"] = 2; options["silent"] = false; options["quad"] = true; OptVal, running_time, status = solve_moment_mip(A, options);

# obj = Dict(); t = Dict(); options = Dict(); options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
# for i = 1:10
#     for car = ["a" "b" "c" "d"]
#         ben = @sprintf("mip_gka%d%s.sparse.mat", i, car)
#         if !isfile(ben)
#             continue
#         end
#         vars = matread(ben);
#         A = vars["L"];
#         for l in [0 4 6 8]
#             options["level"] = l;
#             # OptVal, running_time, status = solve_moment_cert(A, b, c, x00, ϵ, options);
#             OptVal, running_time, status = solve_moment_mip(A, options);
#             obj[@sprintf("%s_level_%d", ben, ql)] = OptVal; t[@sprintf("%s_level_%d", ben, l)] = running_time["solv"];
#         end
#         # OptVal, running_time = solve_moment_cert_tssos(A, b, c, x00, ϵ, options);
#     end
# end

function solve_moment_mip_multi_order(A, options)
    n = size(A, 1); typ = "max";
    L = sparse(diagm(((A.*W)*ones(n,1))[:,1])) - A.*W;
    E0 = sparse([zeros(1+n, n-1) [2 zeros(1,n); ones(n,1) I(n)]]);
    E = sparse([I(n); I(n)]);
    obj = [zeros(1,n+1); [zeros(n,1) A]];
    iter = 0;
    # Shor's relaxation
    NumMom = 1;
    MomConst = Array{Any}(undef, 1);
    MomConst[1] = Dict();
    MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
    NumLoc = 0;
    LocConst = Array{Any}(undef, 0);
    clique = [1:n; 1:n]';
    for i = 1:n
        NumMom = NumMom + 1;
        MomConst = vcat(MomConst, Array{Any}(undef, 1));
        MomConst[NumMom] = Dict();
        MomConst[NumMom]["basis"] = sparse(E[1:n, :]);
        MomConst[NumMom]["ord"] = 0;
        NumLoc = NumLoc + 1;
        LocConst = vcat(LocConst, Array{Any}(undef, 1));
        LocConst[NumLoc] = Dict();
        LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
        LocConst[NumLoc]["basis"] = sparse(I(n));
        LocConst[NumLoc]["typ"] = "==";
        LocConst[NumLoc]["ord"] = 0;
    end
    @printf("\nMAX-CUT Problem: %d vertices, %d edges, order %d, epoch %d\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), options["ord"], iter)
    OptVal, running_time, stat, M = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    while iter <= 10
        iter = iter + 1;
        F = svd(M); U = F.U; u = U[:,1];
        idx = findall(abs.(u.^2-1)==maximum(abs.(u.^2-1)))
        MomConst[idx]["basis"] = sparse(E[clique[idx:idx+iter-1], :]);
        MomConst[idx]["ord"] = options["ord"];
        LocConst[idx]["pol"] = sparse([1 E[idx,:]'; -1 2*E[idx,:]']);
        LocConst[idx]["basis"] = sparse(E[clique[idx:idx+iter-1], :]);
        LocConst[idx]["typ"] = "==";
        LocConst[idx]["ord"] = options["ord"] - 1;
        @printf("\nMAX-CUT Problem: %d vertices, %d edges, order %d, epoch %d\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), options["ord"], iter)
        OptVal, running_time, stat, M = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    end
    @printf("Operation terminated at epoch %d, optimal value is %d, solver status is: ", iter, OptVal); println(stat)
    return OptVal, stat, iter
end
# vars = matread("mip.mat");
# A = vars["L"];
# options = Dict();
# options["ord"] = 2; options["silent"] = true; options["quad"] = true;
# OptVal, running_time, status = solve_moment_mip_multi_order(A, options);

function solve_moment_mip_tssos(A, options)
    d = size(A,1)
    @polyvar x[1:d]
    pop = Array{Any}(undef, 1+d)
    pop[1] = x'*A*x;
    for i = 1:d
        # @printf("%d\n", i)
        pop[i+1] = x[i]*(x[i]-1)
    end
    @printf("\nMixed Integer Programming (MIP): %d variables, %d entries, method %s\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), uppercase(options["method"]))
    running_time = Dict(); OptVal = Dict();
    start = time();
    if options["first"] == true
        OptValFirst, sol, data = blockcpop_first(pop, x, 2, numeq=d, method=options["method"], QUIET=options["silent"]);
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
        n = length(x);
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
                for j = 1:d
                    ssupp[j,i] = MultivariatePolynomials.degree(mon[i],x[j]);
                end
            end
            supp[k] = sparse(ssupp);
        end
        rlorder = ones(Int, m+1);
        rlorder[1] = 2;
        cliques, cql, cliquesize = clique_cdecomp(n, m, supp, rlorder, alg="amd");
        ts = 5;
        mclique,I,ncc,lmc,blocks,cl,blocksize,ub,sizes,ssupp,lt,fbasis,gbasis=get_cblocks_mix(rlorder,m,supp,cliques,cql,cliquesize,ts=ts,method=options["method"],chor_alg="max");
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
# vars = matread("mip.mat");
# A = vars["L"];
# options = Dict();
# options["method"] = "block"; options["silent"] = false; options["first"] = false; options["higher"] = true; options["mix"] = true; options["mix_higher"] = true; OptVal, running_time = solve_moment_mip_tssos(A, options)

function solve_moment_mac(A, options)
    n = size(A, 1); typ = "max";
    E0 = sparse([zeros(1+n, n-1) [2 zeros(1,n); ones(n,1) I(n)]]);
    E = sparse([I(n); I(n)]);
    # L = L/2;
    # for i = 1:n
    #     L[i,i] = L[i,i]/2;
    # end
    # obj = sparse([zeros(n+1,1) E0[1:1+n, n+1:2*n]]);
    # for i = 2:n+1
    #     obj = vcat(obj, [L[i-1:n,i-1] E0[1:n-i+2, n-i+2:2*n-i+1]]);
    # end
    obj = [zeros(1,n+1); [zeros(n,1) A]];
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
            LocConst[i]["typ"] = ">=";
            LocConst[i]["ord"] = 0;
        end
        LocConst = vcat(LocConst, Array{Any}(undef, 1));
        LocConst[n+1] = Dict();
        LocConst[n+1]["pol"] = sparse([-1 zeros(1,n); ones(n,1) E[1:n,:]]);
        LocConst[n+1]["basis"] = sparse(I(n));
        LocConst[n+1]["typ"] = "==";
        LocConst[n+1]["ord"] = 0;
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
                if i+options["level"]-1<=n
                    cl = clique[i:i+options["level"]-1]
                else
                    cl = clique[[1:i+options["level"]-1-n; i:n]]
                end
                NumMom = NumMom + 1;
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[cl, :]);
                MomConst[NumMom]["ord"] = options["ord"];
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[cl, :]);
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
                LocConst[NumLoc]["typ"] = ">=";
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([-1 zeros(1,n); ones(n,1) E[1:n,:]]);
                LocConst[NumLoc]["basis"] = sparse(E[cl, :]);
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
    @printf("\nMAXIMUM-CLIQUE Problem: %d vertices, %d edges, clique %s, order %d, level %d\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), uppercase(options["clique"]), options["ord"], options["level"])
    OptVal, running_time, stat = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    return OptVal, running_time, stat
end
# vars = matread("maxcut_g05_60.0.mat");
# A = vars["A"];
# # n = 10; Random.seed!(1); A = rand(n,n);
# options = Dict();
# options["level"] = 2; options["clique"] = "off"; options["ord"] = 2; options["silent"] = false; options["quad"] = true; OptVal, running_time, status = solve_moment_mac(A, options)

# ["maxcut_g05_60.0.mat" "maxcut_g05_80.0.mat" "maxcut_g05_100.0.mat" "maxcut_pm1d_80.0.mat" "maxcut_pm1d_100.0.mat" "maxcut_pm1s_80.0.mat" "maxcut_pm1s_100.0.mat" "maxcut_pw01_100.0.mat" "maxcut_pw05_100.0.mat" "maxcut_pw09_100.0.mat" "maxcut_w01_100.0.mat" "maxcut_w05_100.0.mat" "maxcut_w09_100.0.mat"]
# obj = Dict(); t = Dict(); options = Dict(); options["range"] = "local"; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
# for ben in ["maxcut_pw05_100.0.mat" "maxcut_pw09_100.0.mat"]
#     vars = matread(ben);
#     A = vars["A"];
#     LB = -Inf; n = size(A,1)
#     for k = 1:1000
#         x = rand(n,1); x = x/sum(x); oo = x'*A*x; LB = maximum([LB oo])
#     end
#     @printf("\nLower bound: %.2f\n", LB)
#     for l in [0 4 6 8]
#         options["level"] = l;
#         OptVal, running_time, status = solve_moment_mac(A, options)
#         obj[@sprintf("%s_level_%d", ben, l)] = OptVal; t[@sprintf("%s_level_%d", ben, l)] = running_time["solv"];
#     end
# end

# obj = Dict(); t = Dict(); options = Dict(); options["range"] = "local"; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
# for n in [60 80 100]
#     for i = 1:2
#         ben = @sprintf("mac_%d_%d.mat", n, i)
#         vars = matread(ben);
#         A = vars["A"];
#         LB = -Inf;
#         for k = 1:1000
#             x = rand(n,1); x = x/sum(x); oo = x'*A*x; LB = maximum([LB oo])
#         end
#         @printf("\nLower bound: %.2f\n", LB)
#         for l in [0 4 6 8]
#             options["level"] = l;
#             OptVal, running_time, status = solve_moment_mac(A, options)
#             obj[@sprintf("%s_level_%d", ben, l)] = OptVal; t[@sprintf("%s_level_%d", ben, l)] = running_time["solv"];
#         end
#     end
# end

function solve_moment_qcqp(A, b, options)
    n = size(A, 1); typ = "min";
    E0 = sparse([zeros(1+n, n-1) [2 zeros(1,n); ones(n,1) I(n)]]);
    E = sparse([I(n); I(n)]);
    # A = 2*A;
    # for i = 1:n
    #     A[i,i] = A[i,i]/2;
    # end
    # obj = sparse([zeros(n+1,1) E0[1:1+n, n+1:2*n]]);
    # for i = 2:n+1
    #     obj = vcat(obj, [A[i-1:n,i-1] E0[1:n-i+2, n-i+2:2*n-i+1]]);
    # end
    obj = [0 b'; [b A]];
    if options["level"] == 0 # Shor's relaxation
        MomConst = Array{Any}(undef, 1);
        MomConst[1] = Dict();
        MomConst[1]["basis"] = sparse(I(n)); MomConst[1]["ord"] = 1;
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
                if i+options["level"]-1<=n
                    cl = clique[i:i+options["level"]-1]
                else
                    cl = clique[[1:i+options["level"]-1-n; i:n]]
                end
                NumMom = NumMom + 1;
                MomConst = vcat(MomConst, Array{Any}(undef, 1));
                MomConst[NumMom] = Dict();
                MomConst[NumMom]["basis"] = sparse(E[cl, :]);
                MomConst[NumMom]["ord"] = options["ord"];
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([1 E[i,:]'; -1 2*E[i,:]']);
                LocConst[NumLoc]["basis"] = sparse(E[cl, :]);
                LocConst[NumLoc]["typ"] = ">=";
                LocConst[NumLoc]["ord"] = options["ord"] - 1;
                NumLoc = NumLoc + 1;
                LocConst = vcat(LocConst, Array{Any}(undef, 1));
                LocConst[NumLoc] = Dict();
                LocConst[NumLoc]["pol"] = sparse([-1 zeros(1,n); ones(n,1) 2*E[1:n,:]]);
                LocConst[NumLoc]["basis"] = sparse(E[cl, :]);
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
    @printf("\nQuadratically Constrained Quadratic Problem (QCQP): %d variables, %d entries, clique %s, order %d, level %d\n", n, (sum(Matrix(1*(A.!=0))) - sum(diag(Matrix(1*(A.!=0)))))/2 + sum(diag(Matrix(1*(A.!=0)))), uppercase(options["clique"]), options["ord"], options["level"])
    OptVal, running_time, stat = solve_moment_manual(typ, obj, MomConst, LocConst, options);
    return OptVal, running_time, stat
end
# vars = matread("mip_bqp50-1.sparse.mat");
# A = vars["L"]; b = zeros(size(A,1),1);
# vars = matread("mip_gka6d.sparse.mat")
# A = vars["L"]; b = zeros(size(A,1),1);
# options = Dict();
# options["level"] = 8; options["clique"] = "off"; options["ord"] = 2; options["silent"] = false; options["quad"] = true; OptVal, running_time, status = solve_moment_qcqp(A, b, options);

# obj = Dict(); t = Dict(); options = Dict(); options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
# # ["mip_bqp50-1.sparse.mat" "mip_bqp100-1.sparse.mat" "mip_gka1a.sparse.mat" "mip_gka2a.sparse.mat" "mip_gka3a.sparse.mat" "mip_gka4a.sparse.mat" "mip_gka5a.sparse.mat" "mip_gka6a.sparse.mat" "mip_gka7a.sparse.mat" "mip_gka8a.sparse.mat" "mip_gka1b.sparse.mat" "mip_gka2b.sparse.mat" "mip_gka3b.sparse.mat" "mip_gka4b.sparse.mat" "mip_gka5b.sparse.mat" "mip_gka6b.sparse.mat" "mip_gka7b.sparse.mat" "mip_gka8b.sparse.mat" "mip_gka9b.sparse.mat" "mip_gka10b.sparse.mat" "mip_gka1c.sparse.mat" "mip_gka2c.sparse.mat" "mip_gka3c.sparse.mat" "mip_gka4c.sparse.mat" "mip_gka5c.sparse.mat" "mip_gka6c.sparse.mat" "mip_gka7c.sparse.mat" "mip_gka1d.sparse.mat" "mip_gka2d.sparse.mat" "mip_gka3d.sparse.mat" "mip_gka4d.sparse.mat" "mip_gka5d.sparse.mat" "mip_gka6d.sparse.mat" "mip_gka7d.sparse.mat" "mip_gka8d.sparse.mat" "mip_gka9d.sparse.mat" "mip_gka10d.sparse.mat" "qcqp_60_1.mat" "qcqp_60_2.mat" "qcqp_80_1.mat" "qcqp_80_2.mat" "qcqp_100_1.mat" "qcqp_100_2.mat"]
# for ben in ["mip_gka9d.sparse.mat" "mip_gka10d.sparse.mat"]
#     vars = matread(ben);
#     A = vars["L"]; b = zeros(size(A,1),1);
#     UB = Inf; n = size(A,1)
#     for k = 1:100000
#         x = rand(n,1); x = x/sqrt(sum(x.^2)); oo = x'*A*x+2*b'x; UB = minimum([UB oo])
#     end
#     @printf("\nUpper bound: %.2f\n", UB)
#     for l in [8]
#         options["level"] = l;
#         OptVal, running_time, status = solve_moment_qcqp(A, b, options)
#         obj[@sprintf("%s_level_%d", ben, l)] = OptVal; t[@sprintf("%s_level_%d", ben, l)] = running_time["solv"];
#     end
# end
# @printf("\n%.4f\n", (539.70+obj["mip_gka9d.sparse.mat_level_8"])/(539.70-130.08))
# @printf("\n%.4f\n", (552.40+obj["mip_gka10d.sparse.mat_level_8"])/(552.40-186.95))

# obj = Dict(); t = Dict(); options = Dict(); options["range"] = "local"; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
# for n in [60 80 100]
#     for i = 1:2
#         ben = @sprintf("qcqp_%d_%d.mat", n, i)
#         vars = matread(ben);
#         A = vars["A"]; b = vars["b"];
#         UB = Inf;
#         for k = 1:1000
#             x = rand(n,1); x = x/sqrt(sum(x.^2)); oo = x'*A*x+2*b'*x; UB = minimum([UB oo])
#         end
#         @printf("\nUpper bound: %.2f\n", UB)
#         for l in [0 4 6 8]
#             options["level"] = l;
#             OptVal, running_time, status = solve_moment_qcqp(A, b, options)
#             obj[@sprintf("%s_level_%d", ben, l)] = OptVal; t[@sprintf("%s_level_%d", ben, l)] = running_time["solv"];
#         end
#     end
# end
