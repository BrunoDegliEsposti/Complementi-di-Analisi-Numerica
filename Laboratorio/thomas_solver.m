function [x] = thomas_solver(c, d, b, e)
%THOMAS Algoritmo di Thomas per la soluzione di sistemi tridiagonali
    n = length(d);
    l = zeros(n,1);
    u = zeros(n,1);
    % c(i) = A(i,i-1) sottodiagonale
    % d(i) = A(i,i)   diagonale
    % b(i) = A(i,i+1) sovradiagonale
    u(1) = d(1);
    for i = 2:n
        l(i) = c(i-1)/u(i-1);
        u(i) = d(i)-l(i)*b(i-1);
    end
    % Sia y=Ux. Allora y = L\e e x = U\y
    x = zeros(n,1);
    y = zeros(n,1);
    y(1) = e(1);
    for i = 2:n
        y(i) = e(i)-y(i-1)*l(i);
    end
    x(n) = y(n)/u(n);
    for i = n-1:-1:1
        x(i) = (y(i)-b(i)*x(i+1)) / u(i);
    end
end








