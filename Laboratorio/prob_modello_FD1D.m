function [u] = prob_modello_FD1D(sigma, f, a, b, ga, gb, N)
%PROB_MODELLO_FD1D Soluzione del problema modello in 1D mediante DF
    h = (b-a)/N;
    x = linspace(a,b,N+1)';
    x_in = x(2:end-1);
    c = -1 * ones(N-2,1) / (h*h);
    d =  2 * ones(N-1,1) / (h*h) + sigma(x_in);
    b = -1 * ones(N-2,1) / (h*h);
    e = f(x_in);
    e(1) = e(1) + ga/(h*h);
    e(end) = e(end) + gb/(h*h);
    u = [ga; thomas_solver(c,d,b,e); gb];
end
