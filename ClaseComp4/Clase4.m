cd('/Users/afl/Dropbox/Universidad/Sexto año/Dynamic Programming/Material_Clases/Códigos/Clase4')

clear
clc

% Generamos la p.d.f. de la distribución Gamma según los parámetros
% descritos:

alpha = 4;
theta = 6;

% Fijamos n:

n = 100;

% Inicializamos un vector para guardar los valores de la p.d.f.

q = zeros(n+1,1);

% Con el loop obtenemos la p.d.f. (podríamos hacerlo de una manera más
% eficiente, no nos vamos a preocupar de ello ahora).

for i = 0:n
    q(i+1,1) = gampdf(i,alpha,theta);
end

% Imponemos los valores de c y de beta:

c = -10;

beta = 0.99;

% El c<0 se interpreta como un "seguro de cesantía", un ingreso al que el
% agente puede acceder en caso de no trabajar.

% Mínimo y máximo del salario:

w_min = 5;
w_max = 20;

% Grilla del salario:

w = linspace(w_min, w_max, n+1);

% Graficamos la distribución de w:

figure(1)
plot(w,q, LineStyle='-', LineWidth=2, Color='red')
title('Densidad de probabilidad de w')
xlabel('w', Interpreter='latex')
ylabel('q(w)', Interpreter='latex')

%% VFI

% Guess inicial

V0 = zeros(n+1,1);

% Número máximo de iteraciones:

itermax = 10000;

% Parámetro de precisión:

epsilon = 1e-7;

for i = 1:itermax

    accion = zeros(n+1,2); % Inicializamos el matriz para guardar los 
                           % valores de las acciones a realizar
    
    % Se guarda el valor de cada posible acción condicional a V_0
    
    for j = 1:n+1
        accion(j,:) = [w(j)/(1-beta), -c + beta * sum(V0.* q)];
    end
    
    % Se encuentra la opción óptima:

    Vn = max(accion, [], 2);

    % Distancia

    d = max(abs(V0 - Vn));

    % Criterio de convergencia:

    if d < epsilon
        V =  Vn;
        break
    else
        V0 = Vn;
    end

end

% El salario de reserva es:

wr = (1-beta) * (- c + beta * sum(V.* q));

% Graficamos

figure(2)
plot(w,V, '-.', LineWidth=2, color='#0c79c4'); hold on
xline(wr, '-', LineWidth=1, color = '#e92427');
title('Función valor')
xlabel('$w$', Interpreter='latex')
ylabel('$v(w)$', Interpreter='latex')

%% Estática comparativa utilizando la función

c_grid = linspace(-15,5,20);

nc = size(c_grid,2);

beta_grid = linspace(0.9,0.999, 20);

nb = size(beta_grid,2);

V_grid = zeros(nb,nc,n+1);
rw_grid = zeros(nb,nc);

for i = 1:nb
    for j = 1:nc
        V = McCallSearch_VFI(beta_grid(i), c_grid(j))';
        V_grid(i,j,:) = V;
        rw_grid(i,j) = (1-beta_grid(i)) * (-c_grid(j) + beta_grid(i) * sum(V' .* q));
    end
end

%%

figure(3)
contourf(c_grid, beta_grid, rw_grid, 'LineStyle', 'none'); % 'LineStyle', 'none' para eliminar las líneas de contorno
xlabel('$c$', Interpreter='latex');
ylabel('$\beta$', Interpreter='latex');
title('Salario de reserva','interpret','latex');
colorbar; % Agregar una barra de color para interpretación visual


figure(4)

surf(c_grid, beta_grid, rw_grid);
colorbar; 
xlabel('c','interpret','latex');
ylabel('$\beta$', Interpreter='latex');
zlabel('Salario de reserva','interpret','latex');
title('Salario de reserva (3D)','interpret','latex');