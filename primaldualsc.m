function [x,y,z,obj,gap,iter,estado] = primaldualsc(c,A,b,x0,y0,z0)
[m,n]=size(A); zer_tol=1.0e-8; tau=0.99995;
if any(b<0) 
    ineg=find(b<0);
    b(ineg)=-b(ineg); A(ineg,:)=-A(ineg,:);
end
c=c(:); b=b(:); A=sparse(A); maxit=150; iter=0; 
if nargin<4 %punto inicial
    xtil=A'*(((A*A')\b)); eps2=100;
    eps1=max([-min(xtil),eps2,norm(b,1)/(eps2*norm(A,1))]);
    x0 = max(xtil,eps1); y0=zeros(m,1); eps3=norm(c,1)+1;
    z0 = max(c,0)+eps3;
end
ind1n=1:n; ind1n=ind1n(:);
x=x0(:); y=y0(:); z=z0(:); estado=0; primal=0; gap=0;
disp('Método de Pontos Interiores Primal-Dual Seguidor de Camino.');
fprintf('Iter\t Primal\t Dual\t gap\t Residuo\n');
while estado == 0
    if iter > maxit
        estado = 3; %alcanza num max de iteraciones
        break;
    end
    iter = iter + 1;
    gam = sum(x.*z);
    if gam < 1
        sig = gam/n; 
    else
        if n < 1000
            sig = 1/n;
        else
            sig = 1/sqrt(n);
        end
    end
    mu =sig*gam/n; d=x./z; d=d(:); D = spconvert([ind1n ind1n d]);
    rp = b - A*x; rd = c - A'*y - z; rc = -x.*z + mu;
    ADAt = A*D*A';
    r1 = rd - rc./x;
    dy = ADAt\(rp + A*(d.*r1));
    dx = d.*(A'*dy - r1);
    dz = (rc - z.*dx)./x;
    ap=1; ad=1;
    bqp = find(dx < -zer_tol); bqd = find(dz < -zer_tol);
    if ~isempty(bqp) 
        ap = tau*min(-x(bqp)./dx(bqp)); 
        ap = min(1,ap); 
    end
    if ~isempty(bqd) 
        ad = tau*min(-z(bqd)./dz(bqd)); 
        ad = min(1,ad); 
    end
    x = x + ap*dx; y = y + ad*dy; z = z + ad*dz;
    primal = sum(c.*x); dual = sum(b.*y);
    gap=max(sum(x.*z),abs(primal-dual));
    gap_rel = gap/(abs(primal)+abs(dual)+1);
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
obj=primal;
