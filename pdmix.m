%Pinturas
A=[6 4 1 0 0;1 2 0 1 0;-1 1 0 0 1]; c=[-5 -4 0 0 0]; b=[24 6 1]; 
lo=zeros(5,1); hi = [1.0e+20 2 1.0e+20 1.0e+20 1.0e+20]; z0=0;
%load('afiro') 
%load('israel')
%%load('maros_r7')
%load('kb2')
%load('fit2p')
[m,n]=size(A);
tic;
if exist('lo','var') ~= 1
    lo=zeros(n,1);
end
if exist('hi','var') ~= 1
    hi = 1.0e+20*ones(n,1);
end
lo = lo(:); hi = hi(:); c=c(:); b=b(:); A = sparse(A);
[cr,Ar,br,lor,hir,row_zeros,col_zeros] = presolving(c,A,b,lo,hi);
Ar = sparse(Ar);
if row_zeros > 0
    fprintf('Presolving: %d filas nulas eliminadas\n',row_zeros);
end
if col_zeros > 0
    fprintf('Presolving: %d filas nulas eliminadas\n',col_zeros);
end
bl = br-Ar*lor; ut = hir - lor;
[ch,Ah,bh,uh,r0,r,s,s0] = rescalplcan(cr,Ar,bl,ut); Ah = sparse(Ah);
tpresc = toc; fprintf('Tiempo de rescalamiento %f segundos\n',tpresc); 
tic;
[xh,vh,yh,zh,wh,objh,gaph,iter,estado] = pdcanmix(Ah,bh,ch,uh);
[x,y,z,v,w,obj,gap] = desrecalplcan(xh,vh,yh,zh,wh,objh,gaph,r0,r,s,s0,uh);
%[x,s,y,z,w,obj,gap,iter,estado] = pdcanmix(A,bl,c,ut);
tpproc = toc;
obj = obj + sum(c.*lo);
x = x + lo;

fprintf('Problema: %d x %d\t nnz=%d\n',m,n,nnz(A));
if row_zeros > 0 | col_zeros > 0
    fprintf('Problema reducido: %d x %d\t nnz=%d\n',m-row_zeros,n-col_zeros,nnz(Ar)); 
end
if estado == 3, fprintf('Alcanzó num max de iteraciones\n'); end
if estado==1
    fprintf('Solución optima obj=%f en %d iteraciones\n',obj,iter);
    fprintf('Tiempo presolving %f y solución %f\n',tpresc,tpproc);
    fprintf('Tiempo total de procesamiento %f segundos\n',tpresc+tpproc);
    if n < 16
        disp(x); disp(obj);
    end
end