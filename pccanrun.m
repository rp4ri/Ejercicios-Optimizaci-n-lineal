clear;
%A=[6 4 1 0 0;1 2 0 1 0;-1 1 0 0 1]; c=[-5 -4 0 0 0]; b=[24 6 1]; 
%lo=zeros(5,1); hi = [1.0e+20 2 1.0e+20 1.0e+20 1.0e+20]; z0=0;
%load('afiro'); 
%load('maros_r7'); 
load('../scsd8c');
[m,n]=size(A); nz=nnz(A); c=c(:); b=b(:); A = sparse(A); rescal = 1;
if exist('lo','var') ~= 1
    lo=zeros(n,1);
end
if exist('hi','var') ~= 1
    hi = 1.0e+20*ones(n,1);
end
lo = lo(:); hi = hi(:); tic;
[c,A,b,lo,hi,nrow_zeros,ncol_zeros] = presolving(c,A,b,lo,hi);
A = sparse(A);
if nrow_zeros > 0
    fprintf('Presolving: %d filas nulas eliminadas\n',nrow_zeros);
end
if ncol_zeros > 0
    fprintf('Presolving: %d filas nulas eliminadas\n',ncol_zeros);
end
b = b-A*lo; ut = hi - lo;
if rescal == 1
    [c,A,b,ut,r0,r,s,s0] = rescalplcan(c,A,b,ut); A = sparse(A);
end
tpresc = toc; fprintf('Tiempo de preprocesamiento %f segundos\n',tpresc);
tic;
[x,v,y,z,w,obj,gap,iter,estado] = pccan(c,A,b,ut);
if rescal == 1
    [x,y,z,v,w,obj,gap] = desrecalplcan(x,v,y,z,w,obj,gap,r0,r,s,s0,ut);
end
tpproc = toc;
obj = obj + sum(c.*lo);
x = x + lo;

fprintf('Problema: %d x %d\t nnz=%d\n',m,n,nz);
if nrow_zeros > 0 | ncol_zeros > 0
    fprintf('Problema reducido: %d x %d\t nnz=%d\n',m-nrow_zeros,n-ncol_zeros,nnz(A)); 
end
if estado == 3, fprintf('Alcanzó num max de iteraciones\n'); end
if estado==1
    fprintf('Solución optima obj=%f en %d iteraciones\n',obj,iter);
    fprintf('Tiempo preprocesamiento %f y solución %f\n',tpresc,tpproc);
    fprintf('Tiempo total de procesamiento %f segundos\n',tpresc+tpproc);
    if n < 16
        disp(x); disp(obj);
    end
end