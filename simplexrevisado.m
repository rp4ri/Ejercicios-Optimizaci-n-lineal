function [x,y,obj,iter,b_idx,n_idx,estado] = simplexrevisado(c,A,b,b_idx,n_idx)
%Resuelve problema PL min c^Tx s.a. Ax=b, x>=0. Ver.2.0 9-09-2020
%Se asume que rango(A)=m=numero de filas de A.
[m,n]=size(A); b=b(:); c=c(:); %vectores columna
%1=optimo, -1=infactble, 2=ilimitado, 0=iternado, -2 base infactible, 3 maxiter
estado = 0; iter = 0; obj=0; maxiter=5000;
y=zeros(m,1); %multiplicadores
x=zeros(n,1); %solucion
zer_tol=10^(-8);
deltab = rand(m,1)*eps; bper=b+deltab;
while estado == 0
   if iter > maxiter
       estado = 3; %por alcazar num max de iteraciones
       break;
   end
   B=A(:,b_idx); N=A(:,n_idx);
   xB=B\bper; %solución básica
   if length(xB(xB<-zer_tol)) > 0
       estado = -2; %base o solucion infactible 
       break;
   end
   iter = iter + 1; %contando iteraciones
   obj = c(b_idx)'*xB; %obj 
   x(b_idx) = xB; x(n_idx)=0;
   y=B'\c(b_idx);
   %if mod(iter,500)==0, fprintf('iteracion=%d obj=%f\n',iter,obj); end
   rN=c(n_idx)-N'*y; %costos reducidos de variables no básicas
   if length(rN(rN<-zer_tol)) == 0 %prueba de optimalidad
       x_opt = B\b; obj=c(b_idx)'*x_opt;
       x(b_idx) = x_opt; x(n_idx)=0;
       estado=1; %optimo
       break; %termina el proceso while
   end
   %Variable de entrada
   [~,q]=min(rN); %q significa la q-ésima entrada de rN
   %variable de salida
   Nq = A(:,n_idx(q));
   Ntq = B\Nq;
   bloqueo = find(Ntq>zer_tol);
   if length(bloqueo) == 0
       estado = 2; %ilimitado
       break; %termina el proceso
   end
   %test de razon
   [~,ind_p] = min(xB(bloqueo)./Ntq(bloqueo));
   p = bloqueo(ind_p);
   aux = b_idx(p); %cambio de base
   b_idx(p) = n_idx(q); %basica p, sale a no básica
   n_idx(q) = aux; %no basica q, entra a la base
end
