function [x,v,y,z,w,obj,gap,iter,estado] = pccan(c,A,b,ut,x,v,y,z,w)
[m,n]=size(A); zer_tol=1.0e-8; eps_tol=1.0e-7; tau=0.99995;
c=c(:); b=b(:); A=sparse(A); u=ut(:); maxit=100; iter=0;
can = ut < 1.0e+20; n1=sum(can);
d=ones(1,n);
ind1n=1:n; ind1n=ind1n(:);
if nargin<5 %punto inicial
    eps2=100; eps3=max(abs(c))+1;
    u = ut(can);
    d(~can)=2; d=d(:); D=spconvert([ind1n ind1n d]);
    ADAt = A*D*A'; AC = sparse(A(:,can));
    if n1 > 0, r=2*b-AC*u; else r=2*b; end 
    g=ADAt\r;
    v = 0.5*(u-AC'*g); x = A'*g; x(can)=x(can) + v; 
    xm=min(x); vm=min(v);
    eps1=max([-min([xm,vm]),eps2,(norm(b,1)+norm(u,1))/(eps2*(norm(A,1)+1))]);
    x = max(x,eps1); v = max(v,eps1);
    z = max(c,0)+eps3; w = z(can); y = zeros(m,1); 
end
estado = 0; obj = 0; x=x(:); v = v(:); y=y(:); z=z(:); w = w(:);
disp('Método de Pontos Interiores Predictor-Corrector.');
fprintf('Iter\t Primal\t Dual\t gap\t Residuo\n');
while estado == 0
    if iter > maxit
        estado = 3; %alcanza num max de iteraciones
        break;
    end
    rp = b - A*x; xz = x.*z; vw = v.*w;
    rd = c-A'*y-z; rd(can) = rd(can) + w; ru = u - v - x(can);
    gam = sum(xz) + sum(vw);
    primal = sum(c.*x); dual = sum(b.*y) - sum(u.*w);
    gap=max(abs(gam),abs(primal - dual));
    gap_rel = gap/(abs(primal) + 1);
    infac_can = norm(ru)/(norm(u)+1);
    infac_primal = norm(rp)/(norm(b)+1);
    infac_dual = norm(rd)/(norm(c)+1);
    max_res = max([gap_rel,infac_primal,infac_can,infac_dual]);
    fprintf('%d\t %0.5f\t %0.5f\t %0.5f\t %0.10f\n',...
        iter,primal,dual,gap,max_res);
    if max_res < eps_tol
        estado = 1; %situación optima
        break;
    end
    iter = iter + 1;
    d = z./x; d(can) = d(can) + w./v; d = 1./d; d=d(:);
    D = spconvert([ind1n ind1n d]);
    ra = -xz; rb = -vw;
    r1 = rd - ra./x; r1(can) = r1(can) + rb./v - (w./v).*ru;
    r = rp + A*(d.*r1);
    ADAt = A*D*A';
    dy = ADAt\r;
    dx = d.*(A'*dy - r1);
    dv = ru - dx(can); dz = (ra - z.*dx)./x; dw = (rb - w.*dv)./v;
    
    ap=1; ad=1;
    if sum(dx<-zer_tol) > 0 
        ap = min(1,tau/max(-dx./x));
    end
    if sum(dv<-zer_tol) > 0 
        ap = min([ap,1,tau/max(-dv./v)]);
    end
    if sum(dz<-zer_tol) > 0
        ad = min(1,tau/max(-dz./z));
    end
    if sum(dw<-zer_tol) > 0
        ad = min([ad,1,tau/max(-dw./w)]);
    end
    gamtil = (x+ap*dx)'*(z+ad*dz) + (v+ap*dv)'*(w+ad*dw);
    sig = (gamtil/gam)^2;
    mu = sig*gam/(n+n1);
    rc = ra + mu - dx.*dz;
    rs = rb + mu - dv.*dw;
    r1 = rd - rc./x; r1(can) = r1(can) + rs./v - (w./v).*ru;
    r = rp + A*(d.*r1);
    dy = ADAt\r;
    dx = d.*(A'*dy - r1);
    dv = ru - dx(can); dz = (rc - z.*dx)./x; dw = (rs - w.*dv)./v;
    
    ap=1; ad=1;
    if sum(dx<-zer_tol) > 0 
        ap = min(1,tau/max(-dx./x));
    end
    if sum(dv<-zer_tol) > 0 
        ap = min([ap,1,tau/max(-dv./v)]);
    end
    if sum(dz<-zer_tol) > 0
        ad = min(1,tau/max(-dz./z));
    end
    if sum(dw<-zer_tol) > 0
        ad = min([ad,1,tau/max(-dw./w)]);
    end
    x = x + ap*dx; v = v + ap*dv; 
    y = y + ad*dy; z = z + ad*dz; w = w + ad*dw;
end
obj=primal;