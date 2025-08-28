clc

%% Priprava (enkrat)
levels = 4;             % št. nivojev za F-cikel
Nfin   = 16;            % št. notranjih točk na finest mreži
omega  = 1;             % SOR-parameter
nu1    = 3; nu2 = 3;    % pred-/post-smooth koraki
maxCycles = 20;         % največ ciklov za primerjavo

% Sestavi Laplace-operator in R/P operaterje za vse nivoje
for k = 1:levels
    n = Nfin/2^(k-1) + 1;
    A{k} = poisson_stencil2D(n);
    if k<levels
        R{k} = restrictionFW2D(n);
        P{k} = interpolation2D((n)/2);
    end
end
A1 = A{1}; R1 = R{1}; P1 = P{1};

% Desna stran (enaka za oba algoritma)
f1 = ones((Nfin)^2,1);

%% 1) Two-grid V-cikel
resTG = zeros(maxCycles,1);
uTG    = zeros(size(f1));
for m = 1:maxCycles
    % V-cikel na finest (k=1)
    % pred-smoothing
    for i=1:nu1
      uTG = relaxGaussSeidel(A1, uTG, f1, omega);
    end
    % restrikcija reziduala
    r = f1 - A1*uTG;               % rezidual
    rc = R1 * r;                   % coarse rezidual
    % coarse solve (en GS korak ali več)
    uc = zeros(size(rc));
    for i=1:nu1+nu2
      uc = relaxGaussSeidel(A{2}, uc, rc, omega);
    end
    % prolongacija in korekcija
    uTG = uTG + P1 * uc;
    % post-smoothing
    for i=1:nu2
      uTG = relaxGaussSeidel(A1, uTG, f1, omega);
    end

    resTG(m) = norm(f1 - A1*uTG);  % shrani L2-normo reziduala
end

%% 2) Four-grid F-cikel
resF = zeros(maxCycles,1);
uF    = zeros(size(f1));
for m = 1:maxCycles
    uF = Fcycle(1, uF, f1, A, R, P, omega, nu1, nu2, levels);
    resF(m) = norm(f1 - A1*uF);
end

%% 3) Primerjava
figure;
semilogy(1:maxCycles, resTG, 'o-', 1:maxCycles, resF, 's-','LineWidth',1.5);
grid on;
xlabel('Število ciklov');
ylabel('‖rezidual‖_2');
legend('Dvomrežni V–cikel','Štirimrežni F–cikel','Location','northeast');
title('Primerjava hitrosti konvergence');