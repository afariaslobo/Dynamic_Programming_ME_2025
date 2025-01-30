% Guía de autoevaluación - Programación Dinámica.
% Otoño 2025.
% Profesor: Eduardo Engel.
% Ayudante: Agustín Farías Lobo.

%% Preliminar

% En esta parte del código definimos el directorio en el que trabajaremos.
% Ello permitirá mantener un código ordenado.

% Primero limpiamos el enviroment ('clear') y la consola ('clc'):

clear
clc

% Para definir el directorio, utilizamos la función 'cd'

cd('/Users/afl/Library/CloudStorage/GoogleDrive-afariaslobo@gmail.com/My Drive/Universidad/Sexto año/Dynamic Programming/Autoevaluación_Matlab')

% Para poder utilizar directamente carpetas dentro del directorio,
% utilizamos la función 'addpath'

addpath('Datos/')
addpath('Funciones/')
addpath('Figuras/')

%% Pregunta 1

% En primer lugar leemos la base de datos con la función readmatrix o
% readtable. Puede ocuparse cualquiera de ellas, pero trabajar con matrices
% puede ser más fácil para ciertas operaciones.

base = readtable("Datos/base_p1.xlsx");
mat = readmatrix("Datos/base_p1.xlsx");

% Primero obtenemos la serie de la inflación con la diferencia del
% logaritmo de IPC con el logaritmo del doceavo rezago del IPC:

pi = log(mat(13:end,3)) - log(mat(1:end-12,3));

% Se crea la función de OLS. Ella está presente en la carpeta "Funciones".

% Definimos la matriz X con las columnas 2 y 4 de la base de datos
% (recordar que perdemos los primeros doce datos para obtener la inflación):

X = mat(13:end,[2,4]);

% Definimos la matriz Y con la inflación:

Y = pi*100;

% Estimamos el modelo con la función OLS.

beta_hat = OLS(X,Y);

% Para obtener las predicciones, basta con multiplicar la matriz X
% (incluyendo la columna de unos) por el vector de coeficientes estimados.

Y_hat = [ones(size(X,1),1), X] * beta_hat;

% Para hacer el gráfico con un formato de fecha, aprovechamos la base en
% formato table. 

tab_plot = table(base.Periodo(13:end), pi*100, Y_hat, ...
                 VariableNames={'Periodo', 'Inflacion', 'Prediccion'});

% A continuación creamos los gráficos:

figure(1)

plot(tab_plot.Periodo, tab_plot.Inflacion, '-', LineWidth=1.5, Color='red'); hold on
% utilizar 'hold on' nos permite sobreponer gráficos
plot(tab_plot.Periodo, tab_plot.Prediccion, '--', LineWidth=1.5, Color='blue');
xlabel('Periodo'); % Ajustamos el nombre del eje x
ylabel('Inflación (%)'); % Ajustamos el nombre del eje y
legend('Inflación observada', 'Inflación predicha', Location='northwest');
% Creamos una leyenda
title('Inflación en doce meses en Chile'); % Incluimos título
subtitle('entre enero de 2014 y noviembre de 2024'); % Incluimos subtítulo

% Guardamos la figura en formato vectorizado (.eps) y en formato
% tradicional (.jpeg).

saveas(1,'Figuras/Figura_P1.eps')
saveas(1,'Figuras/Figura_P1.jpeg')

% Eliminamos las variables en el workspace para comenzar con la siguiente
% pregunta

clearvars

%% Pregunta 2

% En primer lugar establecemos los parámetros fijos:

beta = 0.9;
r = 0.03;

% creamos la grilla para el consumo utilizando linspace:

C_grid = linspace(0,200,20000);

% creamos la grilla de los valores de sigma utilizando linspace:

sigma_grid = linspace(1.1,50,1000);

% Creamos la función U y la dejamos en la carpeta funciones. La llamamos
% U_c.

% Asumimos que C_grid representa c_1. Así, obtenemos c_2 a partir de la
% restricción presupuestaria:

c_2 = (100 - C_grid)*(1+r);

% Para utilizar el loop inicializamos una matriz de ceros:

U = zeros(1000, 20000);

% Ahora usamos el for loop:

for i = 1:1000
    U(i,:) = U_c(C_grid, c_2, sigma_grid(i), beta);
end

% Este loop es sumamente simple: utiliza un índice i que va entre 1 y 100,
% donde en cada step asigna a la fila i el valor de la función de utilidad
% dado C_grid y c_2 (c_1,c_2). Ello es así por como fue creada la función
% U_c, la que se creó de forma que no hubiese problemas al utilizar
% matrices o escalares.

% Ahora utilizamos la función max para encontrar el índice del consumo 
% óptimo dado cada valor de sigma, buscando el valor máximo de la utilidad
% en cada fila (opera sobre las columnas):

[M, C_opt_index] = max(U,[],2); 

% Así, podemos obtener el consumo óptimo para cada periodo con el índice,
% con C_grid y con c_2:

C_opt = [C_grid(C_opt_index); c_2(C_opt_index)];

% Creamos el perfil de consumo:

perfil = C_opt(2,:)./C_opt(1,:);

% Ahora graficamos el perfil de consumo para cada valor de sigma:

figure(2)

plot(sigma_grid, perfil, '-', LineWidth=2, Color='red');
xlabel('$\sigma$', Interpreter='latex'); % Ajustamos el nombre del eje x
                                       % (lo dejamos en formato de .tex
                                       % para utilizar el símbolo de sigma

ylabel('$c_2/c_1$', Interpreter='latex'); % Ajustamos el nombre del eje y
title('Perfil de consumo', Interpreter='latex'); % Incluimos título
subtitle('para diferentes valores de $\sigma$', Interpreter='latex'); % Incluimos subtítulo

% Guardamos la figura en formato vectorizado (.eps) y en formato
% tradicional (.jpeg).

saveas(2,'Figuras/Figura_P2.eps')
saveas(2,'Figuras/Figura_P2.jpeg')

% Eliminamos las variables en el workspace para comenzar con la siguiente
% pregunta

clearvars

%% Pregunta 3

% En primer lugar, la función para calcular derivadas numéricas se define
% en la carpeta de funciones y se nombra der.m.

% La función para generar x_n+1 se encuentra en la carpeta de funciones. Se
% nombra step_n.m.

% La función que genera el algoritmo Newton-Raphson se encuentra en la
% carpeta de funciones. Se nombra NR.m.

% Se encuentran las raíces de las funciones dadas en el enunciado:

f = @(x) x.^2;

g = @(x) x.^3 - (x-2).^5;

h = @(x) exp(-x) + log(x)./log(0.2);

[root_f, iter_f] = NR(f,-10);

[root_g, iter_g] = NR(g,4);

[root_h, iter_h] = NR(h,10);
















