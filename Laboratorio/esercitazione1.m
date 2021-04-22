%% Esercizio 2

c = rand(4,1);
d = rand(5,1);
b = rand(4,1);
A = diag(d,0) + diag(c,-1) + diag(b,1);
e = rand(5,1);
x1 = A\e;
x2 = thomas_solver(c,d,b,e);
norm(x1-x2)

%% Esercizio 3

sigma = @(x) sin(x);
f = @(x) exp(x).*(sin(x).*sin(x)-2*cos(x));
u = @(x) sin(x).*exp(x);
a = 0;
b = pi;
ga = 0;
gb = 0;
N = 50;
u1 = prob_modello_FD1D(sigma,f,a,b,ga,gb,N);
x = linspace(a,b,N+1)';
u2 = u(x);
err = norm(u1-u2,'inf');
% verifica che err tende a zero come h^2 (usa loglog)




