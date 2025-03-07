function [V] = McCallSearch_VFI(beta,c)


% Distribución Gamma para los posibles valores del salario

alpha = 4;
theta = 6;

n = 100;

q = zeros(n+1,1);

for i = 0:n
    q(i+1,1) = gampdf(i,alpha,theta);
end

w_min = 5;
w_max = 20;

w = linspace(w_min, w_max, n+1);

%plot(w,q)

V0 = w'./(1-beta);

itermax = 10000;

epsilon = 1e-7;

for i = 1:itermax

    accion = zeros(n+1,2);

    for j = 1:n+1
        accion(j,:) = [w(j)/(1-beta), -c + beta * sum(V0.* q)];
    end

    Vn = max(accion, [], 2);

    d = max(abs(V0 - Vn));

    if d < epsilon
        V =  Vn;
        break
    else
        V0 = Vn;
    end

end

end