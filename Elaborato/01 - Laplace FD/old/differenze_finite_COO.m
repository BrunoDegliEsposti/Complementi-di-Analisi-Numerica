% x e y sono coordinate cartesiane, rho e theta le relative coordinate
% polari. i e j sono gli indici relativi alla griglia sul piano
% rho-theta che fa da mesh per il metodo delle differenze finite.

Nrho = 50; %max 200, errore 7.78e-7
Ntheta = round(2*pi*Nrho);
f = @(x,y) x./(x.*x+y.*y);
%f = @(x,y) real(1/(x+y*1i));

hrho = (2-1)/Nrho;
htheta = 2*pi/Ntheta;
index = @(i,j) ij_to_index(i,j,Ntheta);
% così non è necessario passare la variabile Ntheta ogni volta

% la griglia ha dimensione (Nrho+1)*Ntheta. A Ntheta manca +1 perché
% la soluzione ha condizione al bordo periodica rispetto a theta
Nvertices = (Nrho+1)*Ntheta;
Nboundary = 2*Ntheta;
Nunknowns = Nvertices - Nboundary;
% i tre vettori che descrivono la matrice sparsa Ah nel formato COO
AhCOO_i = ones(4*Nunknowns, 1);   % indici di riga
AhCOO_j = ones(4*Nunknowns, 1);   % indici di colonna
AhCOO_v = zeros(4*Nunknowns, 1);  % elementi non nulli della matrice
counter = 1;
bh = zeros(Nunknowns, 1);

% iniziamo a riempire Ah con i punti della griglia che, non confinando
% con il bordo, rappresentano il caso generale
for j = 3:Nrho-1
    for i = 1:Ntheta
        % "k" è l'indice di riga in Ah, "l" quello di colonna
        k = index(i,j);
        % il laplaciano in coordinate polari dipende da rho
        rho = 1+(j-1)*hrho;
        % index() gestisce automaticamente la ciclicità di i
        AhCOO_i(counter) = k;
        AhCOO_j(counter) = index(i,j);
        AhCOO_v(counter) = -2/hrho^2 - 2/(rho*htheta)^2;
        counter = counter + 1;
        AhCOO_i(counter) = k;
        AhCOO_j(counter) = index(i-1,j); % theta e i hanno orientazione opposta
        AhCOO_v(counter) = 1/(rho*htheta)^2;
        counter = counter + 1;
        AhCOO_i(counter) = k;
        AhCOO_j(counter) = index(i+1,j);
        AhCOO_v(counter) = 1/(rho*htheta)^2;
        counter = counter + 1;
        AhCOO_i(counter) = k;
        AhCOO_j(counter) = index(i,j-1); % rho e j hanno la stessa orientazione
        AhCOO_v(counter) = 1/hrho^2 - 1/(2*rho*hrho);
        counter = counter + 1;
        AhCOO_i(counter) = k;
        AhCOO_j(counter) = index(i,j+1);
        AhCOO_v(counter) = 1/hrho^2 + 1/(2*rho*hrho);
        counter = counter + 1;
    end
end

% ora occupiamoci dei punti con j=2, che danno contributo anche a bh
for i = 1:Ntheta
        k = index(i,2);
        rho = 1+hrho;
        theta = 2*pi - i*htheta;
        [x_left,y_left] = pol2cart(theta,1);
        dirichlet = f(x_left,y_left);
        l       = index(i,2);
        l_up    = index(i-1,2);
        l_down  = index(i+1,2);
        % l_left non apporta nessun contributo
        l_right = index(i,3);
        Ah(k,l)       = -2/hrho^2 - 2/(rho*htheta)^2;
        Ah(k,l_up)    =  1/(rho*htheta)^2;
        Ah(k,l_down)  =  1/(rho*htheta)^2;
        bh(k)         = -dirichlet * (1/hrho^2 -1/(2*rho*hrho));
        Ah(k,l_right) =  1/hrho^2 + 1/(2*rho*hrho);
end

% passiamo ai punti con j=Nrho
for i = 1:Ntheta
        k = index(i,Nrho);
        rho = 1+(Nrho-1)*hrho;
        theta = 2*pi - i*htheta;
        [x_right,y_right] = pol2cart(theta,2);
        dirichlet = f(x_right,y_right);
        l       = index(i,Nrho);
        l_up    = index(i-1,Nrho);
        l_down  = index(i+1,Nrho);
        l_left  = index(i,Nrho-1);
        % l_right non apporta nessun contributo
        Ah(k,l)       = -2/hrho^2 - 2/(rho*htheta)^2;
        Ah(k,l_up)    =  1/(rho*htheta)^2;
        Ah(k,l_down)  =  1/(rho*htheta)^2;
        Ah(k,l_left)  =  1/hrho^2 - 1/(2*rho*hrho);
        bh(k)         = -dirichlet * (1/hrho^2 +1/(2*rho*hrho));
end

% conversione di Ah dal formato COO a CSC
Ah = sparse(AhCOO_i, AhCOO_j, AhCOO_v, Nunknowns, Nunknowns);
uh = Ah\bh;
u = zeros(size(uh));
X = zeros(size(uh));
Y = zeros(size(uh));
for j = 2:Nrho
    for i = 1:Ntheta
        rho = 1+(j-1)*hrho;
        theta = 2*pi-i*htheta;
        [x,y] = pol2cart(theta,rho);
        u(index(i,j)) = f(x,y);
        X(index(i,j)) = x;
        Y(index(i,j)) = y;
    end
end
e = norm(uh-u, +inf);
scatter3(X,Y,uh-u);
%hold on;
%scatter3(X,Y,uh);







