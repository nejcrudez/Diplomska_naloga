% 
% clear; clc;
% 
% % 1) Ustvarimo Hilbertovo matriko velikosti 3 in desno stran
% n     = 3;
% L     = hilb(n);        % Hilbertova matrika 3×3 (simetrična pozitivno definitna)
% f     = (1:n)';         % npr. f = [1;2;3]
% omega = 1;              % osnovni Gauss–Seidel
% 
% % 2) Začetni približek
% u = zeros(n,1);
% 
% % 3) 5 GS korakov
% fprintf('Iteracija   u (vektor)\n');
% for k = 1:5
%     u = relaxGaussSeidel(L, u, f, omega);
%     fprintf('%2d     [', k);  % izpis na 6 decimalnih mest natančno
%     fprintf(' %8.6f', u);
%     fprintf(' ]\n');
% end
% 
% % 4) Primerjava s točno rešitvijo
% Tocna_resitev = L\f;
% fprintf('Tocna resitev: [');
% fprintf(' %8.6f', Tocna_resitev);
% fprintf(' ]\n');

clear; clc;

% 1) Ustvarimo Hilbertovo matriko velikosti 3 in desno stran
n     = 3;
L     = hilb(n) + eye(n);        % Hilbertova matrika 3×3 (simetrična pozitivno definitna)
f     = (1:n)';         % npr. f = [1;2;3]
omega = 1;              % osnovni Gauss–Seidel

% 2) Začetni približek
u = zeros(n,1);

% 3) 5 GS korakov
fprintf('Iteracija   u (vektor)\n');
for k = 1:5
    u = relaxGaussSeidel(L, u, f, omega);
    fprintf('%2d     [', k);  % izpis na 6 decimalnih mest natančno
    fprintf(' %8.6f', u);
    fprintf(' ]\n');
end

% 4) Primerjava s točno rešitvijo
Tocna_resitev = L\f;
fprintf('Tocna resitev: [');
fprintf(' %8.6f', Tocna_resitev);
fprintf(' ]\n');
