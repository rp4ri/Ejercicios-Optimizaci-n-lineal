function [x,y,z,obj,gap,iter,estado] = pdafinescala(c,A,b,x0,y0,z0)
tau=0.9995; zer_tol=1.0e-5; n = size(A,2);
if any(b<0) %no es estrictamente necesario
    ineg=find(b<0);
    b(ineg)=-b(ineg); A(ineg,:)=-A(ineg,:);
end
if any(x0<0) || any(z0<0)
    fprintf('Su punto inicial no es interior');
    return
end
x=x0(:); y=y0(:); z=z0(:); estado=0; primal=0; gap=0;
c=c(:); b=b(:); A=sparse(A); maxit=100; iter=0;
ind1n=1:n; ind1n = ind1n(:);
disp('Método de Pontos Interiores Primal-Dual Afín Escala.');
fprintf('Iter\t Primal\t Dual\t gap\t Residuo\n');
while estado == 0
    if iter > maxit
        estado=3; %por alcanzar num max de iteraciones
        break; 
    end
    iter = iter + 1;
    rp = b-A*x;
    rd = c-A'*y - z; 
    ra = -x.*z;
    d = x./z; d=d(:);
    D = spconvert([ind1n ind1n d]); %D=Z^{-1}X
    ADAt=A*D*A'; 
    r=rp+A*(d.*(rd-ra./x));
    dy = ADAt\r; %resolucion de ecuaciones normales
    dx = d.*(A'*dy-rd+ra./x);
    dz = (ra-z.*dx)./x;
    ap=1; ad=1; %bloqueadores
    bqx = find(dx < -zer_tol);
    bqz = find(dz < -zer_tol); 
    if ~isempty(bqx) 
        ax = tau*min(-x(bqx)./dx(bqx));
        ap = min(1,ax);
    end
    if ~isempty(bqz) 
        az = tau*min(-z(bqz)./dz(bqz));
        ad = min(1,az);
    end
    apd=min([ap,ad]);
    x = x + apd*dx; 
    y = y + apd*dy; 
    z = z + apd*dz;
    primal = sum(c.*x); 
    dual = sum(b.*y);
    gap = max(sum(x.*z),abs(primal - dual));
    gap_rel = gap/(abs(sum(c.*x)) + abs(sum(b.*y)) + 1);
    infac_primal = norm(b-A*x)/(norm(b)+1);
    infac_dual = norm(c-A'*y-z)/(norm(c)+1);
    max_res = max([gap_rel,infac_primal,infac_dual]);
    fprintf('%d\t %0.5f\t %0.5f\t %0.5f\t %0.10f\n',...
        iter,primal,dual,gap,max_res);
    if max_res < zer_tol
        estado = 1; %situación optima
        break;
    end    
end
obj = primal;