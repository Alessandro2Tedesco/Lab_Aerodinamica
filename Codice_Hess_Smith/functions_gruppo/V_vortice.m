%% Funzione che raccoglie la velocità indotta dalla distribuzione di vortici su ogni pannello in un singolo vettore 

function U_v = V_vortice(N_pannelli, Centro, Estremo_1, Estremo_2, L2G_TransfMatrix, G2L_TransfMatrix, gamma)

U_v = zeros(N_pannelli,2);

for i = 1:N_pannelli

    U_v_old = zeros(1,2);
    Centro_qui = Centro(i, :)';
    
    for j = 1:N_pannelli      
        
        Estremo_1_qui = Estremo_1(j, :)';
        Estremo_2_qui = Estremo_2(j, :)';
    
        L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));      % Prendo la i-esima riga e aggiungo la riga i+1 per ottenere una matrice 2x2
        G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));
    
        U_v_old = U_v_old + ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui)';
    
    end

    U_v(i,:) = U_v_old.*gamma;

end

end
