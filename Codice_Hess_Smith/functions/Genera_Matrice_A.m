function [matriceA] = Genera_Matrice_A(N_pannelli, Centro, Normale, Tangente, Estremo_1, Estremo_2, L2G_TransfMatrix, G2L_TransfMatrix)

matriceA = zeros(N_pannelli+1);

% Calcolo della sottomatrice A_s e del vettore a_v

for i = 1:N_pannelli

    Centro_qui = Centro(i, :)';
    Normale_qui = Normale(i, :)';

        for j = 1:N_pannelli  

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));      % Prendo la i-esima riga e aggiungo la riga i+1 per ottenere una matrice 2x2
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            matriceA(i,j) = dot(ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);        % Prodotto scalare tra il vettore velocità dovuto ad una distribuzione di sorgenti di intensità unitaria su ciascun pannello e versore normale al pannello i-esimo

            matriceA(i,N_pannelli+1) = matriceA(i, N_pannelli+1) + dot(ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);    % Prodotto scalare tra il vettore velocità dovuto ad una distribuzione di vortici di intensità unitaria su ciascun pannello e versore normale al pannello i-esimo 

        end

end


% Calcolo del vettore c_s e di c_v

Centro_Start = Centro(1, :)';
Tangente_Start = Tangente(1, :)';

Centro_End = Centro(end, :)';
Tangente_End = Tangente(end, :)';

c_v = 0;
for j = 1:N_pannelli

    Estremo_1_qui = Estremo_1(j, :)';
    Estremo_2_qui = Estremo_2(j, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(ViSorgente(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);
    c_v = c_v + dot(ViVortice(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);

    a = a + dot(ViSorgente(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);
    c_v = c_v + dot(ViVortice(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);

    matriceA(N_pannelli+1,j) = a;

end

matriceA(N_pannelli+1,N_pannelli+1) = c_v;

end