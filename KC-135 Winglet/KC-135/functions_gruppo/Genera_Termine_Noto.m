function [Termine_Noto] = Genera_Termine_Noto(N_pannelli, Normale, Tangente, U_inf)

% Vettore termine noto sorgenti

Termine_Noto = zeros(N_pannelli+1,1);

for j = 1:N_pannelli
    Normale_qui = Normale(j, :)';
    Termine_Noto(j) = - dot(U_inf, Normale_qui);
end

Tangente_1 = Tangente(1, :)';
Tangente_end = Tangente(end, :)';

% Termine noto vortici

Termine_Noto(N_pannelli+1) = - dot(U_inf, (Tangente_1 + Tangente_end));

end