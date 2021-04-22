f = @(x,y) x./(x.*x+y.*y);
g = @(rho,theta) f(rho*cos(theta),rho*sin(theta));

Nrho = 10; % max 200, errore 6.54e-7
Ntheta = round(2*pi*Nrho); % così hrho è circa htheta
hrho = 1/Nrho;
htheta = 2*pi/Ntheta;
h = sqrt(hrho^2+htheta^2);

N = (Nrho+1)*Ntheta;
tic();
Ah = spalloc(N, N, 5*N); % 5*N perché la molecola 2D ha 5 punti
bh = zeros(N, 1);

% codice relativo ai punti interni alla griglia
for i = 2:Nrho
    for j = 1:Ntheta
        % "k" è l'indice di riga in Ah, "l" quello di colonna
        k = get_k(i,j,Ntheta);
        % il laplaciano in coordinate polari dipende da rho
        rho = 1+(i-1)*hrho;
        % get_k() gestisce automaticamente la ciclicità lungo theta
        l_center = k;
        l_up     = get_k(i,j+1,Ntheta);
        l_down   = get_k(i,j-1,Ntheta);
        l_left   = get_k(i-1,j,Ntheta);
        l_right  = get_k(i+1,j,Ntheta);
        % segue codice poco ottimizzato ma, spero, più chiaro
        Ah(k,l_center) = -2/hrho^2 - 2/(rho*htheta)^2;
        Ah(k,l_up)     =  1/(rho*htheta)^2;
        Ah(k,l_down)   =  1/(rho*htheta)^2;
        Ah(k,l_left)   =  1/hrho^2 - 1/(2*rho*hrho);
        Ah(k,l_right)  =  1/hrho^2 + 1/(2*rho*hrho);
    end
end

% condizioni al bordo di Dirichlet, i = 1
for j = 1:Ntheta
    k = get_k(1,j,Ntheta);
    Ah(k,k) = 1;
    theta = (j-1)*htheta;
    bh(k) = g(1,theta);
end

% condizioni al bordo di Dirichlet, i = Nrho+1
for j = 1:Ntheta
    k = get_k(Nrho+1,j,Ntheta);
    Ah(k,k) = 1;
    theta = (j-1)*htheta;
    bh(k) = g(2,theta);
end

% tempo richiesto per l'assemblaggio
ta = toc();

% tempo richiesto per la soluzione del sistema lineare
tic();
vh = Ah\bh;
ts = toc();

% v contiene la soluzione esatta a precisione di macchina
x = zeros(N,1);
y = zeros(N,1);
v = zeros(N,1);
for i = 1:Nrho+1
    for j = 1:Ntheta
        k = get_k(i,j,Ntheta);
        rho = 1+(i-1)*hrho;
        theta = (j-1)*htheta;
        x(k) = rho*cos(theta);
        y(k) = rho*sin(theta);
        v(k) = g(rho,theta);
    end
end
e = norm(vh-v, +inf);

%%%%% Stima empirica per la stabilità
% norm(inv(full(Ah)),'inf')

%%%%% Plot della soluzione
% fig1=figure(1);
% fig1.Renderer='Painters'; % grafica vettoriale
% scatter3(x,y,vh,12.0,vh,'.')
% axis equal

%%%%% Plot dell'errore
% scatter3(x,y,vh-v,12.0,vh-v,'.')

%%%%% Plot della mesh
% scatter(x,y,8,'black','filled')
% axis([-2.2 2.2 -2.2 2.2])
% axis equal
% grid on








