%% Esercizio 1
% comprendi meglio LDL'
% mesh non uniformi?
% condizioni di neumann?
% scegli un problema di test
% usa mesh sempre piu fini
% fai il disegno come il formaggia

%% Esercizio 2

mu = 1;
gamma = 100;
u = @(x) (exp((gamma/mu)*x)-1)/(exp(gamma/mu)-1);
metodo = input('Scegli un metodo: CFD [0] oppure Upwind [1]:');
switch metodo
    case 0
        mesh = linspace(0,1,1001);
        h = 1e-3;
        Pe = h*gamma/(2*mu);
        lambda = (1+Pe)/(1-Pe);
        u1 = u(mesh);
        u2 = ((lambda^(0:1000))-1)/(lambda^1000-1);
    otherwise
        %
end