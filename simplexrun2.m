%Pinturas
A=[6 4 1 0 0 0;1 2 0 1 0 0;-1 1 0 0 1 0; 0 1 0 0 0 1]; c=[-5 -4 0 0 0 0]; b=[24 6 1 2];
%Problemas de NETLIB
%load('afiro')
%load('scsd8')
[x,y,obj,iter,b_idx,n_idx,estado] = simplex2fases(c,A,b);
