%% Funzione che raccoglie la velocità indotta dalla distribuzione di sorgenti su ogni pannello in un singolo vettore 

function U_s = V_sorgente(N_pannelli, Centro, Estremo_1, Estremo_2, L2G_TransfMatrix, G2L_TransfMatrix, q)

U_s = zeros(N_pannelli,2);
U_s_old = zeros(1,2);

for i = 1:N_pannelli
    
    Centro_qui = Centro(i, :)';

        for j = 1:N_pannelli  

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));      % Prendo la i-esima riga e aggiungo la riga i+1 per ottenere una matrice 2x2
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            U_s_old = U_s_old + ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui)'*q(j);
        
        end

     U_s(i,:) = U_s_old;

end

end