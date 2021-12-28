function [x,y,obj,iter,b_ido,n_ido,estado] = simplex2fases(c,A,b)
[m,n]=size(A); c=c(:); b=b(:); obj=0; x=zeros(n,1); y=zeros(m,1);
iter=0; b_ido=zeros(m,1); n_ido=zeros(n-m,1); estado=0;
bneg=find(b<0);
A(bneg,:) = -A(bneg,:); b=abs(b); %multiplicando por (-1) cuando b_i<0
fprintf('Fase I: ')
AR = [A,eye(m)]; cR=[zeros(1,n) ones(1,m)]; cR=cR(:); 
b_idr=(n+1):(n+m); n_idr=1:n;
tic;
[R,yR,objR,iter1,b_id1,n_id1,estado] = simplexrevisado(cR,AR,b,b_idr,n_idr);
fprintf('Tiempo de procesamiento: %g seg.\n',toc)
if estado ~= 1 
    fprintf('Fase I sin optimo. Terminó con estado=%d\n',estado);
    return
end
if abs(objR) > 0
    fprintf('Problema infactible\n');
    estado = -1;
    return
end
if length(intersect(b_id1,b_idr)) > 0
    fprintf('Hay variables artificiales en la base\n')
    return
end
fprintf('obj=%f en %d iteraciones\n',objR,iter1)
if n < 16
    fprintf('Base óptima:'); disp(b_id1)
end
fprintf('Fase II: ')
n_id2=setdiff(1:n,b_id1);
tic;
[x,y,obj,iter,b_ido,n_ido,estado] = simplexrevisado(c,A,b,b_id1,n_id2);
fprintf('Tiempo de procesamiento: %g seg.\n',toc)
if estado==2
    fprintf('Problema ilimitado\n')
elseif estado == -2
    fprintf('Base infactible\n')
elseif estado == 3
    fprintf('Alcanzó en número máximo de iteraciones\n')
elseif estado==1
    fprintf('obj=%f en %d iteraciones\n',obj,iter)
    if n < 16
        fprintf('Solución Óptima:\n')
        disp(x)
        fprintf('Base Optima:')
        disp(b_ido)
    end
else
    fprintf('Estado desconocido\n')
end