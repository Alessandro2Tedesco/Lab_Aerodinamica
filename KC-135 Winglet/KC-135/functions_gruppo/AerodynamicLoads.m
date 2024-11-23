function [Cl,Cp,Cl_integer,Cm_LE,Cm_c] = AerodynamicLoads(N, U, U_inf, U_inf_normal, gamma, Tangente, Normale, Centro, l, c, x_LE, y_LE)

%% Calcolo del Cl

Circolazione = sum(l)*gamma;        % Circolazione totale sul profilo

% Calcolo del coefficiente di Lift attraverso il Teorema di Jutta-Joukowskilunghezza

Cl = -2*Circolazione/norm(U_inf);           % Il Lift è negativo per via del Th di Kutta-Joukowsky


%% Calcolo del Cp

Cp = zeros(N,1);
for i= 1:N
        Tangente_qui = Tangente(i, :);
        Cp(i) = 1 - dot(U(i,:),Tangente_qui)^2/(norm(U_inf)^2);
end

% Calcolo del coefficiente di Lift attraverso l'integrazione della pressione

Cl_integer = 0;
for i = 1:N
    Normale_qui = Normale(i, :)';
    Cl_integer = Cl_integer + Cp(i)*(l(i)/c)*dot(Normale_qui, U_inf_normal);
end


%% Calcolo del coefficiente di momento rispetto al LE ed 1/4 della corda (valore di Xfoil)

LE_position = [x_LE y_LE 0];
Centro_p = zeros(N,3);
r = zeros(N,3);
for i = 1:N
    Centro_p(i,:) = [Centro(i,1) Centro(i,2) 0];
    r(i,:) = Centro_p(i,:) - LE_position;
end

Normale_p = zeros(N,3);
Cm_LE = 0;
for i = 1:N
    Normale_p(i,:) = [Normale(i,1) Normale(i,2) 0];
    Cm_LE = Cm_LE + Cp(i)*l(i)/(c^2)*dot(cross(r(i,:),Normale_p(i,:)), [0 0 1]);        % Verificare segno sommatoria
end

Cm_c = Cm_LE - Cl/4;      % Xfoil calcola il coefficiente di momento rispetto a 1/4 della corda

end