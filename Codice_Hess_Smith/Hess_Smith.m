%% Hess Smith 

clc
close all
clear all

addpath functions

%% Input

U_inf = 1;                      % Velocità all'infinito [m/s]
alpha = 1;                      % Angolo di incidenza [°]
alpha = 2*pi*alpha/180;         % Angolo di incidenza [rad]

U_inf_x = U_inf * cos(alpha);                           % Componente della velocità asintotica lungo x [m/s]
U_inf_y = U_inf * sin(alpha);                           % Componente della velocità asintotica lungo y [m/s]
U_inf = [U_inf_x; U_inf_y];                             % Vettore velocità asintotica
U_inf_normal = [-U_inf(2); U_inf(1)];                   
U_inf_normal = U_inf_normal ./ norm(U_inf_normal);      % Versore normale alla velocità asintotica

TestCase = 0;

CodiceProfilo = '0012';         % Codice profilo    
Chord = 1;                      % Corda profilo [m]
N_pannelli = 101;                % Numero di pannelli

LE_X_Position = 0;              % Posizione Leading Edge
LE_Y_Position = 0;


%% Creazione profilo

% Per creare il profilo con la discretizzazzione a pannelli vado su Xfoil,
% genero il profilo ed importo i dati

% numero profilo:
% [x,y]=createProfile(CodiceProfilo,N_pannelli,Chord);

Corpo = importXfoilProfile(strcat('NACA_', CodiceProfilo, '.dat'));

% Prima flippa i vettori perché altrimenti se prendessi la prima e la
% seconda colonna non avrei un array, ma una colonna

x = flipud(Corpo.x);
y = flipud(Corpo.y);

Corpo.x = x.*Chord;     % Moltiplico per la corda nel caso questa sia diversa da 1
Corpo.y = y.*Chord;

% Plot profilo discretizzato

figure(1);
plot(Corpo.x, Corpo.y, 'bo-'), grid
xlabel('x')
ylabel('y')
axis equal


%% Creazione di una struttura di pannelli

% Funzione che inserendo la tabella con i dati sul profilo restituisce la
% discretizzazione a pannelli

[Centro, Normale, Tangente, Estremo_1, Estremo_2, beta, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo);


%% Creazione della matrice quadrata A

matriceA = Genera_Matrice_A(N_pannelli, Centro, Normale, Tangente, Estremo_1, Estremo_2, L2G_TransfMatrix, G2L_TransfMatrix);


%% Creazione del termine noto

Termine_Noto = Genera_Termine_Noto(N_pannelli, Normale, Tangente, U_inf);


%% Risoluzione sistema lineare

Soluzione = linsolve(matriceA,Termine_Noto); 

q = Soluzione(1:N_pannelli);
gamma = Soluzione(N_pannelli+1);


%% Calcolo della velocità sui pannelli

% Velocità indotta dalle distribuzioni di sorgenti

U_s = zeros(N_pannelli,2);
U_s_old = zeros(1,2);

for i = 1:N_pannelli

    Centro_qui = Centro(i, :)';
    Normale_qui = Normale(i, :)';

        for j = 1:N_pannelli  

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));      % Prendo la i-esima riga e aggiungo la riga i+1 per ottenere una matrice 2x2
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            U_s_old = U_s_old + ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui)'*q(j);
        end

     U_s(i,:) = U_s_old;

end


% Velocità indotta dalle distribuzioni di vortici

U_v = zeros(N_pannelli,2);
U_v_old = zeros(1,2);

for i = 1:N_pannelli

    Centro_qui = Centro(i, :)';
    Normale_qui = Normale(i, :)';

        for j = 1:N_pannelli  

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));      % Prendo la i-esima riga e aggiungo la riga i+1 per ottenere una matrice 2x2
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            U_v_old = U_v_old + ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui)';
        end

     U_v(i,:) = U_v_old.*gamma;

end


% Campo di velocità

U = zeros(N_pannelli,2);

for i = 1:N_pannelli
    U(i,:) = U_inf' + U_s(i,:) + U_v(i,:);      % Vettore la cui i-esima riga rappresenta il vettore velocità al centro del i-esimo pannello
end


%% Calcolo delle forze aerodinamiche e del Cp

Circolazione = sum(lunghezza)*gamma;        % Circolazione totale sul profilo

% Calcolo del coefficiente di pressione sul i-esimo pannello

Cp = zeros(N_pannelli,1);

for i= 1:N_pannelli
    Tangente_qui = Tangente(i, :)';

    Cp(i) = 1 - (dot(U(i,:),Tangente_qui))^2/(norm(U_inf)^2);
end


% Calcolo del coefficiente di Lift attraverso il Teorema di Jutta-Joukowski

Cl = -2*Circolazione/norm(U_inf);


% Calcolo del coefficiente di Lift attraverso l'integrazione della
% pressione

Cl_pressure_integration = 0;
for i = 1:N_pannelli
    Normale_qui = Normale(i, :)';
    Cl_pressure_integration = Cl_pressure_integration + Cp(i)*(lunghezza(i)/Chord)*dot(Normale_qui, U_inf_normal);
end


% Calcolo del coefficiente di momento rispetto al Leading Edge


