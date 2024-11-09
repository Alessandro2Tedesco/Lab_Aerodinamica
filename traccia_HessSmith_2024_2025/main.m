%% traccia Hess Smith (2024)

clc
close all
clear 

addpath mat_functions

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
NPannelli = 101;                % Numero di pannelli

LE_X_Position = 0;              % Posizione Leading Edge
LE_Y_Position = 0;

%% Creazione profilo

% Per creare il profilo con la discretizzazzione a pannelli vado su Xfoil,
% genero il profilo ed importo i dati

% numero profilo:
% [x,y]=createProfile(CodiceProfilo,NPannelli,Chord);

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
        
%% Inizializzazione matrici e vettori

% Ora che ho i pannelli, posso inizializzare la matrice ed i vettori

NCols = NPannelli + 1;
NRows = NCols;
matriceA = zeros(NRows, NCols);
TermineNoto = zeros(NRows, 1);

%% Creazione della matrice quadrata As e del vettore a_v


for i = 1:NPannelli
    index_i = i;            % Riga

    Centro_qui = Centro(i, :)';
    Normale_qui = Normale(i, :)';

    indexStart_colonna = 0;

        for j = 1:NPannelli
            index_j = indexStart_colonna + j;   % Colonna

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));      % Prendo la i-esima riga e aggiungo la riga i+1 per ottenere una matrice 2x2
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            matriceA(index_i, index_j) = dot(ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);        % Prodotto scalare tra il vettore velocità dovuto ad una distribuzione di sorgenti di intensità unitaria su ciascun pannello e versore normale al pannello i-esimo

            matriceA(index_i, NPannelli+1) = matriceA(index_i, NPannelli+1) + dot(ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);    % Prodotto scalare tra il vettore velocità dovuto ad una distribuzione di vortici di intensità unitaria su ciascun pannello e versore normale al pannello i-esimo 

        end

end


%% Creazione delle componenti dei vettori c_s e c_v

% Trovo il centro ed il versore tangente del primo e dell'ultimo pannello

Centro_Start = Centro(1, :)';
Tangente_Start = Tangente(1, :)';

Centro_End = Centro(end, :)';
Tangente_End = Tangente(end, :)';

b = 0;
for j = 1:NPannelli
    
    index_j = j;

    Estremo_1_qui = Estremo_1(j, :)';
    Estremo_2_qui = Estremo_2(j, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(ViSorgente(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);
    b = b + dot(ViVortice(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);

    a = a + dot(ViSorgente(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);
    b = b + dot(ViVortice(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);

    matriceA(NPannelli+1, index_j) = a;

end

matriceA(NPannelli+1, NPannelli+1) = b;



%% Creazione del termine noto

% Vettore termine noto sorgenti

TermineNoto = zeros(NPannelli+1,1);

for j = 1:NPannelli
    index = j;
    
    Normale_qui = Normale(j, :)';

    TermineNoto(index) = - dot(U_inf, Normale_qui);
end

Tangente_1 = Tangente(1, :)';
Tangente_end = Tangente(end, :)';

% Termine noto vortici

TermineNoto(NPannelli+1) = - dot(U_inf, (Tangente_1 + Tangente_end));


%% Risoluzione sistema lineare

Soluzione = linsolve(matriceA,TermineNoto); 

q = Soluzione(1:NPannelli);
gamma = Soluzione(NPannelli+1);


%% Calcolo del cp e della velocità sui pannelli

%%


