function [u] = nonlinear_CFD1D(F, Fy, Fz, a, b, ga, gb, N)
    % F,Fy,Fz: R3 -> R
    h = (b-a)/N;
    x = linspace(a,b,N+1)';
    % applicazione del metodo di newton per un sistema N-1 x N-1
    y = linspace(ga, gb, N+1)'; %approssimiamo la soluzione con una retta
    % definizione di phi, di cui cerchiamo uno zero
    phiy = zeros(N+1,1);
    for i = 2:N
        phiy(i) = (-y(i-1)+2*y(i)-y(i+1))/(h*h) + F(x(i),y(i),(y(i+1)-y(i-1))/(2*h));
    end
    % i vettori Jl, Jd e Ju contengono le tre diagonali della matrice
    % jacobiana di phi (calcolata in y)
    Jl = zeros(N-2,1);
    Jd = zeros(N-1,1);
    Ju = zeros(N-2,1);
    successo = 0;
    for k = 1:100
        for i = 2:N
            Jl(i) = -1/(h*h)-(1/(2*h))*Fz(x(i),y(i),;
            Jd(i-1) = 0;
            Ju(i) = 0;
        end
        Jd(N-1) = 0;
        dy = [0; thomas_solver(Jl, Jd, Ju, phiy(2:end-1)); 0];
        y = y - dy;
        for i = 2:N
            phiy(i) = (-y(i-1)+2*y(i)-y(i+1))/(h*h) + F(x(i),y(i),(y(i+1)-y(i-1))/(2*h));
        end
        if norm(dy,'inf') < 1e-10 * (1+norm(y,'inf'))
            successo = 1;
            break
        end
    end
    if successo == 0
        error('Il metodo di newton ha fallito');
    end
    u = y;
end

