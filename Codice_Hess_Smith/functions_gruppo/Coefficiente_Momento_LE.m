function Cm_LE = Coefficiente_Momento_LE(N_pannelli, Cp, l, c, Centro, x_LE)

r = zeros(N_pannelli,1);
for i = 1:N_pannelli
    r(i) = Centro(i,1) - x_LE;
end

Cm_LE = 0;
for i = 1:N_pannelli
    Cm_LE = -Cp(i)*l(i)/(c^2)*r(i);
end

end