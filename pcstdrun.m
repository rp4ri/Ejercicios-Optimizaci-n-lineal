%Pinturas
%A=[6 4 1 0 0 0;1 2 0 1 0 0;-1 1 0 0 1 0;0 1 0 0 0 1]; c=[-5 -4 0 0 0 0]; b=[24 6 1 2]; 
%load('afiro'); 
%load('maros_r7'); 
%load('../scsd8'); 
load('osa_60'); %load('nug07')
%load('D:/repo/optim/LP_MATLAB/osa_60')
[m,n]=size(A);
c=c(:); b=b(:); A = sparse(A);
tic;
[x,y,z,obj,gap,iter,estado] = pcstd(c,A,b);
tpproc = toc;

fprintf('Problema: %d x %d\t nnz=%d\n',m,n,nnz(A));
if estado == 3, fprintf('Alcanzó num max de iteraciones\n'); end
if estado==1
    fprintf('Solución optima obj=%f en %d iteraciones\n',obj,iter);
    fprintf('Tiempo total de procesamiento %f segundos\n',tpproc);
    if n < 16
        disp(x); disp(obj);
    end
end