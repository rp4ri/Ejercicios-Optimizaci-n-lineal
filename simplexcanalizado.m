function [x,y,obj,iter,b_idx,n_idx,estado] = simplexcanalizado(c,A,b,lb,ub,b_idx,n_idx)
b = b(:); c=c(:); lb=lb(:); ub=ub(:);
global deltab;
[m,n] = size(A); A=sparse(A);
estado = 0; iter = 0; obj=0; maxiter=5000;
y=zeros(m,1); %multiplicadores
x=zeros(n,1); %solucion
zer_tol=1.0e-5;
if isempty(deltab) || length(deltab)~= m 
 deltab = rand(m,1)*eps;
end
bper=b+deltab;
x(n_idx(n_idx>0)) = lb(n_idx(n_idx>0)); 
x(-n_idx(n_idx<0)) = ub(-n_idx(n_idx<0));
while estado == 0
   if iter > maxiter
       estado = 3; %por alcazar num max de iteraciones
       break;
   end
   B=sparse(A(:,b_idx)); N=sparse(A(:,abs(n_idx)));
   bder = bper-N*x(abs(n_idx));
   xB = B\bder; %sol basica
   ilbB = find(lb(b_idx) > -1.0e20);
   iubB = find(ub(b_idx) < 1.0e20);

   if any(xB(ilbB)<lb(b_idx(ilbB))-zer_tol) || any(xB(iubB)>ub(b_idx(iubB))+zer_tol)
       estado = -2; %base o solucion infactible 
       break;
   end 
   iter = iter + 1; %contando iteraciones
   x(b_idx) = xB;
   obj = c'*x; %obj 
   y=B'\c(b_idx);
   %if mod(iter,500)==0, fprintf('iteracion=%d obj=%f\n',iter,obj); end
   rN=sign(n_idx(:)).*(c(abs(n_idx))-N'*y); %costos reducidos de variables no básicas
   %disp(rN);
   if length(rN(rN<-zer_tol)) == 0 %prueba de optimalidad
       bder = b-N*x(abs(n_idx));
       x_opt = B\bder; x(b_idx) = x_opt; obj=c'*x;
       estado=1; %optimo
       break; %termina el proceso while
   end
   %Variable de entrada
   [~,q]=min(rN); %q significa la q-ésima entrada de rN
   %variable de salida
   Nq = sign(n_idx(q))*A(:,abs(n_idx(q)));
   Ntq = B\Nq;
   maxeps = abs(ub(abs(n_idx(q)))-lb(abs(n_idx(q)))); cambio=0; %cambio de lim
   bqlb = intersect(find(Ntq >= zer_tol),ilbB);
   if ~isempty(bqlb)
       [maxep2,ind_p] = min((xB(bqlb)-lb(b_idx(bqlb)))./Ntq(bqlb));
       if maxep2 < maxeps
           maxeps = maxep2; cambio=1; %sale al limite inferior
           p = bqlb(ind_p); salea=lb(b_idx(p));
       end
   end
   bqub = intersect(find(Ntq<-zer_tol),iubB);
   if ~isempty(bqub)
       [maxep3,ind_p] = min((xB(bqub)-ub(b_idx(bqub)))./Ntq(bqub));
       if maxep3 < maxeps
           maxeps = maxep3; cambio=-1; %sale al limite superior
           p = bqub(ind_p); salea=ub(b_idx(p));
       end
   end
   if maxeps >= 1.0e20
       estado = 2; %ilimitado
       break; %termina el proceso
   end
   if abs(cambio) > 0 %cambio de base
       x(b_idx(p)) = salea;
       aux = b_idx(p); %cambio de base
       b_idx(p) = abs(n_idx(q)); %basica p, sale a no básica
       n_idx(q) = sign(cambio)*aux; %no basica q, entra a la base
   else %cambio de límite
       if n_idx(q)<0 
           x(-n_idx(q)) = lb(-n_idx(q));
       else
           x(n_idx(q)) = ub(n_idx(q));
       end
       n_idx(q) = -n_idx(q); %cambio de límite
   end
end
