function [x,v,y,z,w,obj,gap,iter,estado] = pdcanmix(A,b,c,ut,x0,v0,y0,z0,w0)
[m,n]=size(A); zer_tol=1.0e-8; eps_tol=1.0e-7; tau=0.99995;
c=c(:); b=b(:); A=sparse(A); ut=ut(:); maxit=200; iter=0;
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
    vtil = 0.5*(u-AC'*g); xtil = A'*g; xtil(can)=xtil(can) + vtil; 
    xm=min(xtil); vm=min(vtil);
    eps1=max([-min([xm,vm]),eps2,(norm(b,1)+norm(u,1))/(eps2*(norm(A,1)+1))]);
    x0 = max(xtil,eps1); v0=max(vtil,eps1);
    z0 = max(c,0)+eps3; w0=z0(can); y0=zeros(m,1); 
end
estado = 0; obj = 0;
x=x0(:); v=v0(:); y=y0(:); z=z0(:); w=w0(:);
disp('Método de Pontos Interiores Primal-Dual Canalizado Seguidor de Camino Mixto.');
fprintf('Iter\t Primal\t Dual\t gap\t Residuo\n');
while estado == 0
    if iter > maxit
        estado = 3; %alcanza num max de iteraciones
        break;
    end
    rp = b - A*x; xtz = x.*z; vtw = v.*w;
    ru = u - x(can) - v;
    rd = c-A'*y-z; rd(can)=rd(can)+w;
    gam = sum(xtz) + sum(vtw);
    primal = sum(c.*x); dual = sum(b.*y)-sum(u.*w);
    gap=max(abs(gam),abs(primal - dual));
    gap_rel = gap/(abs(primal) + 1);
    infac_primal = norm(rp)/(norm(b)+1);
    infac_dual = norm(rd)/(norm(c)+1);
    infac_can = norm(ru)/(norm(u)+1);
    max_res = max([gap_rel,infac_primal,infac_dual,infac_can]);
    fprintf('%d\t %0.5f\t %0.5f\t %0.5f\t %0.10f\n',...
        iter,primal,dual,gap,max_res);
    if max_res < eps_tol
        estado = 1; %situación optima
        break;
    end
    iter = iter + 1;
    if gam < 1
        sig = gam/(n+n1); 
    else
        if n < 1000
            sig = 1/(n+n1);
        else
            sig = 1/sqrt(n+n1);
        end
    end
    %sig = 1/sqrt(n+n1);
    mu = sig*gam/(n+n1);
    d = z./x; d(can) = d(can)+ w./v; d = 1./d; d=d(:);
    D = spconvert([ind1n ind1n d]);
    rc = -xtz + mu; rs = -vtw + mu;
    r1 = rd - rc./x; r1(can) = r1(can) + rs./v - (w./v).*ru;
    r = rp + A*(d.*r1);
    ADAt = A*D*A';
    dy = ADAt\r;
    dx = d.*(A'*dy - r1);
    dv = ru - dx(can); dz = (rc - z.*dx)./x; dw = (rs - w.*dv)./v;
    
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
obj=primal;