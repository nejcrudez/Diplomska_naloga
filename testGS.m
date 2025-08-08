% testGS_Hilbert5.m
clear; clc;

% 1) Sestavi Hilbertovo matriko in desno stran
n     = 3;
L     = hilb(n);        % Hilbertova matrika 5×5 (SPD)
f     = (1:n)';         % npr. f = [1;2;3;4;5]
omega = 1;              % Gauss–Seidel

% 2) Začetni približek
u = zeros(n,1);

% 3) Naredi nekaj GS korakov
fprintf('Iteracija   u (vektor)\n');
for k = 1:5
    u = relaxGaussSeidel(L, u, f, omega);
    fprintf('%2d     [', k);
    fprintf(' %8.6f', u);
    fprintf(' ]\n');
end

% 4) Primerjava s točno rešitvijo
Tocna_resitev = L\f;
fprintf('Tocna resitev: [');
fprintf(' %8.6f', Tocna_resitev);
fprintf(' ]\n');
