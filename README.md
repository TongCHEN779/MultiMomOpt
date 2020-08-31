# MultiMomOpt
Multi-order and sub-level hierarchy based on standard Lasserre's hierarchy

## Dependencies
JuMP, AMD, LinearAlgebra, MAT, SparseArrays, LightGraphs, GraphPlot, DynamicPolynomials, MosekTools, Printf, Dualization, Random

## Usage
Maximize x1*x2, subject to x1^2+x2^2<=1.

### Do it manually
```Julia
using MultiMomOpt
typ = "max"; obj = [1 1 1];
options = Dict();
options["silent"] = true; options["quad"] = false;
MomConst = Array{Any}(undef, 1);
MomConst[1] = Dict();
MomConst[1]["basis"] = [1 0; 0 1];
MomConst[1]["ord"] = 2;
LocConst = Array{Any}(undef, 1);
LocConst[1] = Dict();
LocConst[1]["basis"] = [1 0; 0 1];
LocConst[1]["pol"] = [1 0 0; -1 2 0; -1 0 2];
LocConst[1]["typ"] = ">=";
LocConst[1]["ord"] = 1;
OptVal, running_time, status = solve_moment_manual(typ, obj, MomConst, LocConst, options);
```

### Do it automatically
```Julia
using MultiMomOpt
@polyvar x[1:2];
obj = x[1]*x[2];
var = x; typ = "max";
MomConst = Array{Any}(undef, 1);
MomConst[1] = Dict(); MomConst[1]["basis"] = x; MomConst[1]["ord"] = 2;
LocConst = Array{Any}(undef, 1);
LocConst[1] = Dict(); LocConst[1]["basis"] = x; LocConst[1]["pol"] = 1-sum(x.^2); LocConst[1]["typ"] = ">="; LocConst[1]["ord"] = 1;
options = Dict(); options["silent"] = true; options["quad"] = false;
OptVal, running_time, status = solve_moment_auto(typ, var, obj, MomConst, LocConst, options);
```

## Examples
### MAXCUT problem
#### Statement
#### Usage
```Julia
using MultiMomOpt
vars = matread("maxcut.mat");
A = vars["A"]; W = ones(size(A, 1), size(A, 1));
options = Dict();
options["level"] = 3; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
OptVal, running_time, status = solve_moment_maxcut(A, W, options)
```

### MAX-CLIQUE problem
#### Statement
#### Usage
```Julia
using MultiMomOpt
vars = matread("mac_5.mat");
A = vars["A"];
options = Dict();
options["level"] = 3; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
OptVal, running_time, status = solve_moment_mac(A, options)
```

### Mixed Integer Programming (MIP)
#### Statement
#### Usage
```Julia
using MultiMomOpt
vars = matread("mip.mat");
A = vars["L"];
options = Dict();
options["level"] = 15; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
OptVal, running_time, status = solve_moment_mip(A, options);
```

### QCQP problem
#### Statement
#### Usage
```Julia
using MultiMomOpt
vars = matread("qcqp_5_1.mat");
A = vars["A"]; b = vars["b"];
options = Dict();
options["level"] = 3; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
OptVal, running_time, status = solve_moment_qcqp(A, b, options)
```

### Lipschitz Constant Estimation problem (one hidden layer)
#### Statement
#### Usage
```Julia
using MultiMomOpt
vars = matread("lip_test.mat");
A = vars["A"]; b = vars["b"]; c = vars["c"]; x00 = vars["x00"]; eps = 0.1;
options = Dict();
options["range"] = "global"; options["level"] = 6; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
OptVal, running_time, status = solve_moment_lip_one_layer(A, b, c, x00, eps, options);
```

### Robustness Certification problem
#### Statement
#### Usage
```Julia
using MultiMomOpt
vars = matread("lip_test.mat");
A = vars["A"]; b = vars["b"]; c = vars["c"]; x00 = vars["x00"]; eps = 0.1;
options = Dict();
options["range"] = "global"; options["level"] = 0; options["clique"] = "off"; options["ord"] = 2; options["silent"] = true; options["quad"] = true;
OptVal, running_time, status = solve_moment_cert(A, b, c, x00, eps, options);
```

## References
The Lipschitz Constant Estimation problem is referred to [Semialgebraic Optimization for Lipschitz Constants of ReLU Networks](https://arxiv.org/abs/2002.03657). The Robustness Certification problem is referred to [Semidefinite relaxations for certifying robustness to adversarial examples](https://arxiv.org/abs/1811.01057). For more information, contact me: tchen@laas.fr.

# Ellipsoid Problem
Finding the minimum volume ellipsoid of the image under ReLU networks

## Dependencies
Plots, SumOfSquares, DynamicPolynomials, MosekTools, LinearAlgebra, Printf, MathOptInterface, StatsBase

## Usage
```Julia
n = 10; Q10 = Matrix(-I(n)); b10 = zeros(n,1); c10 = 1;
L10 = 2; P10 = Array{Any}(undef, L10); ξ10 = Array{Any}(undef, L10);
for i = 1:L10-1
    P10[i] = randn(n,n)*0.5; ξ10[i] = randn(n,1);
    for j = 1:n
        P10[i][j,j] = 1
    end
end
P10[L10] = rand(2,n)*0.5; ξ10[L10] = rand(2,1); P10[L10][1,1] = 1; P10[L10][2,2] = 1;

p11, m = OuterApproximationPlotSampling(Q10, b10, c10, 1, P10, ξ10, "Morari", L10);
p21 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=2, morari=m, meth="cycle_v");
p31 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=2, morari=m, meth="order_v");
p12 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=3, morari=m, meth="order_v");
p22 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=4, morari=m, meth="order_v");
p32 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=5, morari=m, meth="order_v");
p13 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=6, morari=m, meth="order_v");
p23 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=7, morari=m, meth="order_v");
p33 = OuterApproximationPlotSampling(Q10, b10, c10, 2, P10, ξ10, "sublevel", L10, lv=8, morari=m, meth="order_v");

plot(p11, p21, p31, p12, p22, p32, p13, p23, p33, layout = Plots.grid(3,3), fmt = :png)
```

## References
[Reach-SDP: Reachability Analysis of Closed-Loop Systems with Neural Network Controllers via Semidefinite Programming](https://arxiv.org/abs/2004.07876). For more information, contact me: tchen@laas.fr.
