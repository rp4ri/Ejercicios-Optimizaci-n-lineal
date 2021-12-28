function [x,v,y,z,w,obj,gap,iter,estado] = pdsccan(A,bcan,c,lo,hi,x0,v0,y0,z0,w0)
[m,n]=size(A); zer_tol=1.0e-6; eps_tol=1.0e-8; tau=0.99995;
if isempty(lo)
    lo=zeros(n,1);
end
if isempty(hi)
    hi = 1.0e+20*ones(n,1);
end
lo = lo(:); hi = hi(:); bcan=bcan(:);
b = bcan-A*lo; u = hi - lo;
if any(b<0) %no es estrictamente necesario
    ineg=find(b<0);
    b(ineg)=-b(ineg); A(ineg,:)=-A(ineg,:);
end
c=c(:); b=b(:); A=sparse(A); u=u(:); maxit=150; iter=0; 
if nargin<6 %Estrategia de puntos iniciales
    AAt = A*A';
    g=AAt\(2*b-A*u); eps2=100;
    vtil = 0.5*(u-A'*g); xtil = 0.5*(u+A'*g); xm=min(xtil); vm=min(vtil);
    eps1=max([-min(xm,vm),eps2,(norm(b,1)+norm(u,1))/(eps2*(norm(A,1)+1))]);
    x0 = max(xtil,eps1); v0=max(vtil,eps1); y0=zeros(m,1); 
    eps3=norm(c,1)+1; z0 = max(c,0)+eps3; w0=z0;
end
ind1n=1:n; ind1n=ind1n(:); estado = 0; obj = 0;
x=x0(:); v=v0(:); y=y0(:); z=z0(:); w=w0(:);
disp('Método de Puntos Interiores Primal-Dual Canalizado Seguidor de Camino.');
fprintf('Iter\t Primal\t Dual\t gap\t Residuo\n');
while estado == 0
    if iter > maxit
        estado = 3; %alcanza num max de iteraciones
        break;
    end
    rp = b - A*x; ru = u - x - v; rd = c - A'*y - z + w; 
    ra = -x.*z; rb = -v.*w;
    gam = sum(-ra) + sum(-rb);
    primal = sum(c.*x); dual = sum(b.*y)-sum(u.*w);
    gap=max(abs(gam),abs(primal - dual));
    gap_rel = gap/(abs(primal) + 1);
    infac_primal = norm(rp)/(norm(b)+1);
    infac_dual = norm(rd)/(norm(c)+1);
    infac_can = norm(ru)/(norm(u)+1);
    max_res = max([gap_rel,infac_primal,infac_dual,infac_can]);
    fprintf('%d\t%0.5f\t%0.5f\t%0.5f\t%0.10f\n',...
        iter,primal,dual,gap,max_res);
    if max_res < eps_tol
        estado = 1; %situación optima
        break;
    end
    iter = iter + 1;
    %gam = sum(x.*z) + sum(v.*w);
    if gam < 1
        sig = gam/(2*n); 
    else
        if n < 1000
            sig = 1/(2*n);
        else
            sig = 1/sqrt(2*n);
        end
    end
    %sig = 1/sqrt(2*n);
    mu =sig*gam/(2*n); d=1./(z./x+w./v); d=d(:); D = spconvert([ind1n ind1n d]);
    %rp = b - A*x; ru = u - x - v; rd = c - A'*y - z + w; 
    rc = ra + mu; rs = rb + mu;
    ADAt = A*D*A';
    r1 = rd - rc./x + rs./v - (w./v).*ru;
    r = rp + A*(d.*r1);
    dy = ADAt\r;
    dx = d.*(A'*dy - r1);
    dv = ru - dx; dz = (rc - z.*dx)./x; dw = (rs - w.*dv)./v;
    
    ap=1; ad=1;
    bqx = find(dx < -zer_tol); bqv = find(dv < -zer_tol);
    bqz = find(dz < -zer_tol); bqw = find(dw < -zer_tol);
    if ~isempty(bqx) 
        ax = tau*min(-x(bqx)./dx(bqx));
        ap = min(1,ax);
    end
    if ~isempty(bqv) 
        av = tau*min(-v(bqv)./dv(bqv)); 
        ap = min(ap,av); 
    end
    if ~isempty(bqz) 
        az = tau*min(-z(bqz)./dz(bqz));
        ad = min(1,az);
    end
    if ~isempty(bqw) 
        aw = tau*min(-w(bqw)./dw(bqw));
        ad = min(ad,aw);
    end
    x = x + ap*dx; v = v + ap*dv; y = y + ad*dy; 
    z = z + ad*dz; w = w + ad*dw;
end
obj = primal + sum(c.*lo);
x = x + lo;