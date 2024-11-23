function [Centro, Normale, Tangente, Estremo_1, Estremo_2, beta, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo_Input)

NPannelli = length(Corpo_Input.x)-1;

Centro = zeros(NPannelli, 2);                   % Vettore che contiene le coordinate x,y del centro di ciascun pannello
Normale = zeros(NPannelli, 2);                  % Vettore che contiene sulle righe il versore normale del i-esimo pannello
Tangente = zeros(NPannelli, 2);                 % Vettore che contiene sulle righe il versore tangente del i-esimo pannello
Estremo_1 = zeros(NPannelli, 2);                % Vettore che contiene le coordinate x,y del primo estremo di ogni pannello
Estremo_2 = zeros(NPannelli, 2);                % Vettore che contiene le coordinate x,y del secondo estremo di ogni pannello
beta = zeros(NPannelli, 1);                     % Inclinazione di ciascun pannello rispetto al sistema di riferimento globale
L2G_TransfMatrix = zeros(NPannelli, 2, 2);      % Matrice di rotazione che consente di passare dal sistema di riferimento locale a quello globale
G2L_TransfMatrix = zeros(NPannelli, 2, 2);      % Matrice di rotazione che consente di passare dal sistema di riferimento globale a quello locale
lunghezza = zeros(NPannelli,1);                 % Lunghezza in modulo del i-esimo pannello

for i = 1:NPannelli
    
    Centro(i, 1) = (Corpo_Input.x(i) + Corpo_Input.x(i+1))/2;
    Centro(i, 2) = (Corpo_Input.y(i) + Corpo_Input.y(i+1))/2;
    
    Estremo_1(i, 1) = Corpo_Input.x(i);
    Estremo_1(i, 2) = Corpo_Input.y(i);
    
    Estremo_2(i, 1) = Corpo_Input.x(i+1);
    Estremo_2(i, 2) = Corpo_Input.y(i+1);
    
    dy = Estremo_2(i, 2) - Estremo_1(i, 2);
    dx = Estremo_2(i, 1) - Estremo_1(i, 1);
    
    angle = atan2(dy, dx);
    
    if (abs(angle) < 1e-12)
        angle=0; 
    end
    
    beta(i) = angle;
    sinAngle = sin(angle);
    cosAngle = cos(angle);
    
    if(abs(sinAngle) < 1e-12)
        sinAngle = 0;
    end
    
    if(abs(cosAngle) < 1e-12)
        cosAngle = 0;
    end
    
    L2G_TransfMatrix(i, :, :) = [cosAngle ,  -sinAngle;
                                 sinAngle,  cosAngle];
                             
    G2L_TransfMatrix(i, :, :) = [cosAngle ,  sinAngle;
                                 -sinAngle,  cosAngle];
                             
    Normale(i, 1) = -sinAngle;
    Normale(i, 2) = cosAngle;
    
    Tangente(i, 1) = cosAngle;
    Tangente(i, 2) = sinAngle;
    
    lunghezza(i) = norm(Estremo_2(i, :) - Estremo_1(i, :));
    
end

end

