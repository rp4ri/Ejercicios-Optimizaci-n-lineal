A=[6 4 1 0 0;1 2 0 1 0;-1 1 0 0 1]; c=[-5 -4 0 0 0]; b=[24 6 1]; lb = zeros(6,1); ub = [Inf;2;Inf;Inf;Inf];
%load('scsd8c'); lb=zeros(2750,1); ub=u;
%load('afiro');lb=lo; ub=hi;
%load('kb2'); lb=lo; ub = hi;
[x,y,obj,iter,b_idx,n_idx,estado] = simplexcan2fases(c,A,b,lb,ub);
