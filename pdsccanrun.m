%Pinturas
A=[6 4 1 0 0;1 2 0 1 0;-1 1 0 0 1]; c=[-5 -4 0 0 0]; b=[24 6 1]; lo=zeros(5,1); 
hi = [1.0e+20 2 1.0e+20 1.0e+20 1.0e+20]; z0=0;
%load('pinturascan')
%load('afiro') 
%load('israel')
%load('maros_r7')
%load('kb2') %canalizado
%load('fit2p')
%load('../can5x2')
%A=[1 -2 -3 4;-2 4 -4 -5]; b=[-10 10]; c=[2 3 -5 -4];
[m,n]=size(A);
tic;
[x,s,y,z,w,obj,gap,iter,estado] = pdsccan(A,b,c,lo,hi);
tpproc = toc;
fprintf('%d x %d\t nnz=%d\t z0=%f\n',m,n,nnz(A),z0);
if estado == 3, fprintf('Alcanzó num max de iteraciones\n'); end
if estado==1
    fprintf('Solución optima obj=%f en %d iteraciones\n',obj,iter);
    fprintf('Tiempo de procesamiento %f segundos\n',tpproc);
    if n < 16
        disp(x); disp(obj);
    end
end