function [x,y,z,obj,gap,iter,estado] = pcstd(c,A,b,x,y,z)
[m,n]=size(A); zer_tol=1.0e-8; eps_tol=1.0e-8; tau=0.99995;
c=c(:); b=b(:); A=sparse(A); maxit=100; iter=0;
ind1n=1:n; ind1n=ind1n(:);
if nargin<4 %punto inicial
    eps2=100; eps3=max(abs(c))+1;
    AAt = A*A'; x=A'*(AAt\b);
    eps1=max([-min(x),eps2,(norm(b,1))/(eps2*(norm(A,1)))]);
    x = max(x,eps1); z = max(c,0)+eps3; y = zeros(m,1); 
end
estado = 0; obj = 0; x=x(:); y=y(:); z=z(:);
disp('Método de Puntos Interiores Predictor-Corrector.');
fprintf('Iter\t Primal\t Dual\t gap\t Residuo\n');
while estado == 0
    if iter > maxit
        estado = 3; %alcanza num max de iteraciones
        break;
    end
    rp = b - A*x; xz = x.*z;
    rd = c-A'*y-z;
    gam = sum(xz);
    primal = sum(c.*x); dual = sum(b.*y);
    gap=max(abs(gam),abs(primal - dual));
    gap_rel = gap/(abs(primal) + 1);
    infac_primal = norm(rp)/(norm(b)+1);
    infac_dual = norm(rd)/(norm(c)+1);
    max_res = max([gap_rel,infac_primal,infac_dual]);
    fprintf('%d\t %0.5f\t %0.5f\t %0.5f\t %0.10f\n',...
        iter,primal,dual,gap,max_res);
    if max_res < eps_tol
        estado = 1; %situación optima
        break;
    end
    iter = iter + 1;
    d = x./z; d=d(:);
    D = spconvert([ind1n ind1n d]);
    ra = -xz;
    r1 = rd - ra./x;
    r = rp + A*(d.*r1);
    ADAt = A*D*A';
    dy = ADAt\r;
    dx = d.*(A'*dy - r1);
    dz = (ra - z.*dx)./x;
    
    ap=1; ad=1;
    if sum(dx<-zer_tol) > 0 
        ap = min(1,tau/max(-dx./x));
    end
    if sum(dz<-zer_tol) > 0
        ad = min(1,tau/max(-dz./z));
    end
    gamtil = (x+ap*dx)'*(z+ad*dz);
    sig = (gamtil/gam)^2;
    mu = sig*gam/n;
    rc = -xz + mu - dx.*dz;
    r1 = rd - rc./x;
    r = rp + A*(d.*r1);
    dy = ADAt\r;
    dx = d.*(A'*dy - r1);
    dz = (rc - z.*dx)./x;
    
    ap=1; ad=1;
    if sum(dx<-zer_tol) > 0 
        ap = min(1,tau/max(-dx./x));
    end
    if sum(dz<-zer_tol) > 0
        ad = min(1,tau/max(-dz./z));
    end
    x = x + ap*dx; y = y + ad*dy; z = z + ad*dz;
end
obj=primal;