%Calculamos el ?rea y el per?metro de un rectangulo de lados 
% a y b.
% A = area = a x b
% P = Perimetro = 2(a+b)
% Autor: Porfirio Su?agua S. 13 de septiembre 2021
% Materia MAT-258
function [A,P] = AreaPerimetro(a,b)
A = 0; P=0;
if a<0
    error('la medida del lado a debe ser positivo');
end
if b<0, error('el lado b debe ser positivo'); end
%C?lulo del area
A = productonum(a,b);
%C?lcuo del per?metro
P = 2*(a+b);
fprintf('El area del rectangulo es A=%f\n',A);
fprintf('El perimetro del rectangulo es P=%0.8f\n',P);
end

function y = productonum(a,b)
y=a*b;
end