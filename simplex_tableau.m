%c = [-5; -4; 0; 0; 0; 0];
%A = [6 4 1 0 0 0 ; 1 2 0 1 0 0 ; -1 1 0 0 1 0 ; 0 1 0 0 0 1];
%b = [24; 6; 1; 2];
%A=[1 3 0 1 1 0 0; 2 1 0 0 0 1 0; 1 0 4 1 0 0 1];
%c=[2; 4; 1; 1; 0; 0; 0]; b=[4; 3; 3];
%A=[5 -5 13 -2 1; 1 -2 5 -1 1]; c=[1; 6; -7; 1; 5];b=[20; 8];
c=[4;1;0;0]; b=[3;6;4]; A=[3 1 0 0;4 3 -1 0;1 2 0 1];
[opt_x, tobj]= simplex_method(c, A, b);
disp(opt_x); disp(tobj);

function isnoneg=isnonnegative(x)
    isnoneg = isempty( x( x < 0) );
end

function [b_idx, x_B, B]=initial_BFS(A, b)
    [m, n] = size(A);

    comb = nchoosek(1:n, m);
    for i=size(comb,1):-1:1
      b_idx = comb(i,:);
      B = A(:, b_idx);
      x_B = B\b;
      if isnonnegative(x_B)
        return;
      end
    end
    error('Infactible')
end

function print_tableau(t)
    [m, n] = size(t.Y);

    hline0 = repmat('-', 1,6);
    hline1 = repmat('-', 1,7*n);
    hline2 = repmat('-', 1,7);
    hline = [hline0, '+', hline1, '+', hline2];

    disp(hline)

    fprintf('%6s|', '');
    for j=1:length(t.z_c)
      fprintf('%6.2f ', t.z_c(j));
    end
    fprintf('| %6.2f\n', t.obj);

    disp(hline)

    for i=1:m
      fprintf('x[%2d] |', t.b_idx(i));
      for j=1:n
        fprintf('%6.2f ', t.Y(i,j));
      end
      fprintf('| %6.2f\n', t.x_B(i));
    end

    disp(hline)
end

function t=pivoting(t)
    [m, n] = size(t.Y);
    [entering, exiting] = pivot_point(t);
    fprintf('Pivoteo: entrada = x_%d, salida = x_%d\n',entering,t.b_idx(exiting));

    % Pivoting: exiting-row, entering-column
    % updating exiting-row
    coef = t.Y(exiting, entering); 
    t.Y(exiting, :) =t.Y(exiting, :)./coef;
    t.x_B(exiting) =t.x_B(exiting)./coef;

    % updating other rows of Y
    for i = setdiff(1:m, exiting)
      coef = t.Y(i, entering);
      t.Y(i, :) = t.Y(i, :) - coef * t.Y(exiting, :);
      t.x_B(i) = t.x_B(i) - coef * t.x_B(exiting);
    end

    % updating the row for the reduced costs
    coef = t.z_c(entering); 
    t.z_c = t.z_c - coef * t.Y(exiting, :);
    t.obj = t.obj - coef * t.x_B(exiting);

    % Updating b_idx 
    t.b_idx(exiting) = entering;
    %OK t.b_idx(find(t.b_idx==t.b_idx(exiting))) = entering;
end

function [entering, exiting] = pivot_point(t)
    % Finding the entering variable index
%     enter = find(t.z_c < 0,1); %Regla de Bland
%     if isempty(enter)
%       error("Optimo")
%     end

    [mincost, enter] = min(t.z_c); %Regla de Dantzig
    if mincost >= 0
      error('Optimo')
    end

    entering = enter(1);%enter(2);
    % min ratio test / finding the exiting variable index
    pos_idx = find( t.Y(:, entering) > 0 );

    if isempty(pos_idx)
      error('Ilimitado');
    end
    %disp(pos_idx); disp(entering);
    [~,imin] = min( t.x_B(pos_idx) ./ t.Y(pos_idx, entering) );
    exiting = pos_idx(imin);
end

function t=initialize(c, A, b)
    t = struct();
    [m, n] = size(A);

    % Finding an initial BFS
    [b_idx, x_B, B] = initial_BFS(A,b);

    Y = B\A;
    c_B = c(b_idx);
    %obj = -dot(c_B, x_B)
    obj = -c_B'*x_B;

    % z_c is a row vector
    z_c = zeros(1,n);
    n_idx = setdiff(1:n, b_idx);
    %z_c[n_idx] = c_B' * inv(B) * A[:,n_idx] - c[n_idx]'
    z_c(n_idx) = c(n_idx)' - c_B' * (B\A(:,n_idx));
    t.z_c = z_c;
    t.Y = Y;
    t.x_B = x_B;
    t.obj = obj;
    t.b_idx = b_idx; 
end

function costoopt = isOptimal(t)
    costoopt = isempty(find( t.z_c < 0,1 ));
end

function [opt_x, tobj] = simplex_method(c, A, b)
    t = initialize(c, A, b);
    print_tableau(t);

    while ~isOptimal(t)
      t=pivoting(t);
      print_tableau(t);
    end

    opt_x = zeros(length(c),1);
    opt_x(t.b_idx) = t.x_B;
    tobj = -t.obj;
end