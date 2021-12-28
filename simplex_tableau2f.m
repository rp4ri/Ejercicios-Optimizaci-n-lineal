% c = [-5; -4; 0; 0; 0; 0];b = [24; 6; 1; 2];
% A = [6 4 1 0 0 0 ; 1 2 0 1 0 0 ; -1 1 0 0 1 0 ; 0 1 0 0 0 1];

c=[4;1;0;0]; b=[3;6;4]; A=[3 1 0 0;4 3 -1 0;1 2 0 1];

% A=[1 3 0 1 1 0 0; 2 1 0 0 0 1 0; 1 0 4 1 0 0 1];
% c=[2; 4; 1; 1; 0; 0; 0]; b=[4; 3; 3];

%A=[5 -4 13 -2 1; 1 -1 5 -1 1]; c=[1; 6; -7; 1; 5];b=[20; 8];

if length(b(b < 0)) > 0
    for i=1:m
        if b(i) < 0
            b(i) = -b(i);
            A(i,:) = -A(i,:);
        end
    end
end
[opt_x,tobj,nits] = metodo_simplex(c, A, b);
disp(opt_x); disp(tobj); disp(nits);
 
  function ispos=ispositive(x)
    ispos = ~isempty( x( x > 1.0e-4) );
  end

  function [b_one, nit] = phase_one(A,b)
    [m, n] = size(A);
    cR = zeros(n+m,1);
    basis = (n+1):(n+m);
    cR(basis)=ones(m,1);
    AR = sparse([A eye(m,m)]);

    tableau = initialize(cR, AR, b, basis);
    print_tableau(tableau);

    [tableau,b_one,nit] = simplex_iterations(tableau);
    x_opt = zeros(length(cR),1);
    x_opt(tableau.b_idx) = tableau.x_B;
    x_r = x_opt(basis);

    if ispositive(x_r)
        error('Problema Infactible')
    else
        if length(intersect(basis,b_one)) > 0
            error('Variables artificiales en la base')
        else
            return;
        end
    end
  end

  function print_tableau(t)
    formatorac=1; %0 para real
    [m, n] = size(t.Y);
    if n > 15
        return
    end
    hline0 = repmat('-', 1,6);
    if formatorac
        hline1 = repmat('-', 1,6*n);
    else
        hline1 = repmat('-', 1,7*n);
    end
    hline2 = repmat('-', 1,7);
    hline = [hline0, '+', hline1, '+', hline2];

    disp(hline);

    fprintf('%6s|', '');
    for j=1:length(t.z_c)
        if formatorac
            fprintf('%6s',strtrim(rats(t.z_c(j))));
        else
            fprintf('%6.2f ', t.z_c(j));
        end
    end
    if formatorac
        fprintf('| %6s\n',strtrim(rats(t.obj)));
    else
        fprintf('| %6.2f\n', t.obj);
    end

    disp(hline);

    for i=1:m
      fprintf('x[%2d] |', t.b_idx(i))
      for j=1:n
          if formatorac
            fprintf('%6s',strtrim(rats(full(t.Y(i,j)))));
          else
            fprintf('%6.2f ', full(t.Y(i,j)));
          end
      end
      if formatorac
        fprintf('| %6s\n',strtrim(rats(t.x_B(i))));
      else
        fprintf('| %6.2f\n', t.x_B(i));
      end
    end

    disp(hline);
  end

  function t=pivoting(t)
    [m, n] = size(t.Y);
    [entering, exiting] = pivot_point(t);
    t.nit = t.nit + 1;
    if n <= 15
        fprintf('Pivoteo: entrada = x_%d, salida = x_%d\n',entering,t.b_idx(exiting));
    end
    % Pivoting: exiting-row, entering-column
    % updating exiting-row
    coef = t.Y(exiting, entering);
    t.Y(exiting, :) = t.Y(exiting, :)/coef;
    t.x_B(exiting) = t.x_B(exiting)/coef;

    % updating other rows of Y
    for i=setdiff(1:m, exiting)
      coef = t.Y(i, entering);
      t.Y(i, :) = t.Y(i, :)- coef * t.Y(exiting, :);
      t.x_B(i) = t.x_B(i) - coef * t.x_B(exiting);
    end

    % updating the row for the reduced costs
    coef = t.z_c(entering);
    t.z_c = t.z_c - coef * t.Y(exiting, :);
    t.obj = t.obj - coef * t.x_B(exiting);

    % Updating b_idx 
    t.b_idx(exiting) = entering;
    %t.b_idx[findall(t.b_idx.==t.b_idx[exiting])] .= entering
  end

  function [entering, exiting]=pivot_point(t)
    % Finding the entering variable index
    % Enter criteria Bland=1, Dantzig=2
    method = 2;
    if method == 1
        enter = find(t.z_c < -1.0e-4,1); %Regla de Bland
        if isempty(enter)
          error('Optimo');
        end
    elseif method == 2
        [mincost, enter] = min(t.z_c); %Regla de Dantzig
        if mincost >= -1.0e-4
          error('Optimo');
        end
    end
    
    entering = enter(1);

    % min ratio test / finding the exiting variable index
    pos_idx = find( t.Y(:, entering) > 1.0e-4 );
    
    if isempty(pos_idx)
      error('Ilimitado');
    end
    %exiting = pos_idx[ argmin( t.x_B[pos_idx] ./ t.Y[pos_idx, entering] ) ]
    [~,imin] = min( t.x_B(pos_idx) ./ t.Y(pos_idx, entering) );
    exiting = pos_idx(imin);
  end

  function t=initialize(c, A, b, basis)
    t = struct();
    [m, n] = size(A);

    b_idx = basis;
    B = sparse(A(:, b_idx));
    x_B = B\b;
    
    Y = B\A;
    A = sparse(A);
    Y = sparse(Y);
    c_B = c(b_idx);
    %obj = -dot(c_B, x_B)
    obj = -c_B'*x_B;

    % z_c is a row vector
    z_c = zeros(1,n);
    n_idx = setdiff(1:n, b_idx);
    z_c(n_idx) = c(n_idx)' - c_B' * Y(:,n_idx); 
    fprintf('Problema de %d x %d\n',m,n);
    t.z_c = z_c;
    t.Y = Y;
    t.x_B = x_B;
    t.obj = obj;
    t.b_idx = b_idx; 
    t.nit = 0;
  end

  function isopt=isOptimal(tableau)
    isopt = isempty(find( tableau.z_c < -1.0e-4,1 ));
  end

  function [tableau,base,nit] = simplex_iterations(tableau)

    while ~isOptimal(tableau)
      tableau=pivoting(tableau);
      print_tableau(tableau)
    end
    base = tableau.b_idx;
    nit = tableau.nit;
  end

  function [opt_x,tobj,nits]=metodo_simplex(c, A, b)
    disp('Phase I');
    tic;
        [b_one, nit] = phase_one(A,b);
    t1=toc;
    fprintf('N?mero de iteraciones Fase I: %2d\n', nit)
    fprintf('Tiempo de procesamiento Fase I: %f\n',t1);

    disp('Phase II');
    tic;
        tableau = initialize(c, A, b, b_one);
        print_tableau(tableau);
        [tableau,b_two,nits] = simplex_iterations(tableau);
    t2=toc;
    fprintf('Tiempo de procesamiento Fase II: %f\n',t2);
    
    opt_x = zeros(length(c),1);
    opt_x(tableau.b_idx) = tableau.x_B;
    fprintf('N?mero de iteraciones Fase II: %2d\n', tableau.nit);
    if length(b_two)<15, fprintf('Base optima:'); disp(b_two); end
    tobj = -tableau.obj;
  end