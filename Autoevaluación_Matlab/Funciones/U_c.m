function [U] = U_c(c1,c2, sigma, beta)

% recordar que para aplicar una operación componente por componente de una
% matriz, basta con incluir un . antes del operador.

U = c1.^(1-sigma)/(1-sigma) + beta * c2.^(1-sigma)/(1-sigma); 

U(c1 <= 0) = -10^(10);
U(c2 <= 0) = -10^(10);

end