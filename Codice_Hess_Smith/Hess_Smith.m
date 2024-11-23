%% Hess Smith 

clc
close all
clear 

addpath functions
addpath functions_gruppo

%% Input

U_inf = 1;                      % Velocità all'infinito [m/s]
alpha = 2;                      % Angolo di incidenza [°]
alpha = pi*alpha/180;           % Angolo di incidenza [rad]

U_inf_x = U_inf * cos(alpha);                           % Componente della velocità asintotica lungo x [m/s]
U_inf_y = U_inf * sin(alpha);                           % Componente della velocità asintotica lungo y [m/s]
U_inf = [U_inf_x; U_inf_y];                             % Vettore velocità asintotica
U_inf_normal = [-U_inf(2); U_inf(1)];                   
U_inf_normal = U_inf_normal ./ norm(U_inf_normal);      % Versore normale alla velocità asintotica

TestCase = 0;

CodiceProfilo = '0012';         % Codice profilo    
Chord = 1;                      % Corda profilo [m]
N_pannelli = 101;               % Numero di pannelli

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

figure(1)
plot(Corpo.x, Corpo.y, 'bo-'), grid
title("Profilo")
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

U_s = V_sorgente(N_pannelli, Centro, Estremo_1, Estremo_2, L2G_TransfMatrix, G2L_TransfMatrix, q);


% Velocità indotta dalle distribuzioni di vortici

U_v = V_vortice(N_pannelli, Centro, Estremo_1, Estremo_2, L2G_TransfMatrix, G2L_TransfMatrix, gamma);


% Campo di velocità

U = zeros(N_pannelli,2);
for i = 1:N_pannelli
    U(i,:) = U_inf' + U_s(i,:) + U_v(i,:);      % Vettore la cui i-esima riga rappresenta il vettore velocità al centro del i-esimo pannello
end


%% Calcolo dei coefficienti aerodinamici e del Cp

[Cl,Cp,Cl_integer,Cm_LE,Cm_c] = AerodynamicLoads(N_pannelli, U, U_inf, U_inf_normal, gamma, Tangente, Normale, Centro, lunghezza, Chord, LE_X_Position, LE_Y_Position);

fprintf('Il Coefficiente di Lift del profilo è pari a: %f \n', Cl)

% Plot coefficiente di pressione

Cp_ventre = -Cp(1:round(N_pannelli/2));
Cp_dorso = -Cp(round(N_pannelli/2):end);
x_ventre = x(1:round(N_pannelli/2));
x_dorso = flip(x(1:round(N_pannelli/2)));

figure(2)
plot(x_ventre, Cp_ventre, 'ro-'), grid
hold on
plot(x_dorso, Cp_dorso, 'bo-')
hold off
title("Distribuzione Cp")
legend("Cp_v_e_n_t_r_e", "Cp_d_o_r_s_o")
xlabel("x")
ylabel("-Cp")

fprintf("Il coefficiente di momento rispetto al Leading Edge è pari a: %f \n", Cm_LE)


%% Calcolo angolo di Theodorsen