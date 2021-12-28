function [x,y,obj,iter,b_idx,n_idx,estado] = simplexcan2fases(c,A,b,lb,ub)
global deltab;
if nargin < 5 
    ub = Inf*ones(length(c),1);
    if nargin < 4 
        lb = zeros(length(c),1); 
    end;
end
if isempty(lb) || all(lb==0)
    lb = zeros(length(c),1);
end
if isempty(ub)
    ub = Inf*ones(length(c),1);
end;
b=b(:); c=c(:); lb=lb(:); ub=ub(:);
if any(b<0) %cambia de signos a las restricciones con lado derecho negativo
    ineg=find(b<0); 
    b(ineg)=-b(ineg); A(ineg,:)=-A(ineg,:);
end
[m,n] = size(A); %zer_tol = 1.0e-5; 
A = sparse(A); v = zeros(n,1); 
x=v; y=zeros(m,1); iter=0; obj=0; b_idx=[]; n_idx=[];
disp('Fase-I');
%identificando los limites inferior, superior e variables libres.
indlb = find(lb>-1.0e20); v(indlb) = lb(indlb); libre = setdiff(1:n,indlb);
indub = libre(find(ub(libre)<1.0e20)); v(indub) = ub(indub); 
if isempty(indub) 
    n_idr = indlb; 
else
  n_idr = [indlb -indub]; libre = setdiff(libre,indub);
end;
numlib = length(libre); %existencia de variables libres.
if length(libre) > 0
    disp('No se admiten variables libres');
    return;
else %todos tienen cotas
  b_idr = (n+1):(n+m);
  AR = [A sparse(1:m,1:m,sign(b-A*v+eps*ones(length(b),1)))];
end;

%definición de limitantes para las variables "artificiales"
lbr = [lb; zeros(m,1)]; ubr = [ub; Inf*ones(m,1)]; cR = [zeros(n,1); ones(m,1)];
fprintf('Modelo PL canalizado con %d variables, %d artificiales y %d restricciones\n',n,m,m);
if size(b_idr,1)>1, b_idr=b_idr'; end
if size(n_idr,1)>1, n_idr=n_idr'; end
tic;
[xR,y,objR,iter1,b_id1,n_id1,estado] = simplexcanalizado(cR,AR,b,lbr,ubr,b_idr,n_idr);
tp1=toc; %para calcular o tempo de procesamentoç
fprintf('Tiempo de procesamiento Fase-I: %g seg.\n',toc)
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
    fprintf('Hay variables artificiales en la base\n');
    disp('Se podría hacer un proceso extra para proceguir, pero aun no implementado');
    return
end
fprintf('obj=%f en %d iteraciones\n',objR,iter1)
if n < 16
    fprintf('Base óptima:'); disp(b_id1); 
    fprintf('No Básicas:'); disp(n_id1);
end

if objR > 0 %cR'*xR > zer_tol
    disp('El problema PL es infactible'); 
    return
end
disp('Fase-II');
x=zeros(n,1); y=zeros(m,1); iter=0; obj=0;
fprintf('Modelo PL canalizado con n = %d variables y m = %d restricciones\n',n,m);
if(~isempty(libre)), fprintf('Variables libres: %d',length(libre)); return; end
n_id2= n_id1(abs(n_id1) <= n);
tic; %para calcular o tempo de procesamento da Fase-II
[x,y,obj,iter,b_idx,n_idx,estado] = simplexcanalizado(c,A,b,lb,ub,b_id1,n_id2);
tp2=toc;
fprintf('Tiempo de procesamiento Fase-II: %g seg.\n',tp2)
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
        fprintf('Base Optima:'); disp(b_idx); 
        fprintf('No Básicas:'); disp(n_idx);
    end
else
    fprintf('Estado desconocido\n')
end
