%% multigrid_Fcycle_4grid.m
clear; clc;

%---- PARAMETRI ----%
levels = 4;         % število nivojev (4-mrežni)
Nfin = 16;          % št. notranjih točk na finest mreži (mora biti deljivo s 2^(levels-1))
omega = 1;          % SOR parameter
nu1   = 3;          % pred-glajenje
nu2   = 3;          % po-glajenje

%---- PRIPRAVA OPERATORJEV ----%
% A{k} … Laplace na nivoju k
% R{k} … restrikcija iz nivoja k -> k+1
% P{k} … interpolacija iz nivoja k+1 -> k
for k = 1:levels
    n = Nfin/2^(k-1) + 1;       % število točk na stran na nivoju k
    A{k} = poisson_stencil2D(n);
    if k < levels
        R{k} = restrictionFW2D(n);
        P{k} = interpolation2D((n-1)/2);
    end
end

%---- DESNA STRAN IN ZAČETNI PRIBLIŽEK ----%
f1 = ones(Nfin^2,1);    % npr. f≡1 na finest
u0 = zeros(size(f1));

%---- EN KRATKI F-CIKEL ----%
u = Fcycle(1, u0, f1, A, R, P, omega, nu1, nu2, levels);

%---- PRIKAZ REŠITVE ----%
U = reshape(u, Nfin, Nfin);
figure
surf(U);
title('4-mrežni F-cikel – končna rešitev');
xlabel('i'); ylabel('j'); zlabel('u_{i,j}');


%%---------------------------------------------------------------------------%%
function u = Fcycle(level, u, f, A, R, P, omega, nu1, nu2, maxLevel)
% Rekurzivni F-cikel
% level     … trenutni nivo (1 = finest)
% u, f      … vektor rešitve in desna stran na tem nivoju
% A, R, P   … celice z operatorji
% omega     … SOR parameter
% nu1, nu2  … št. pre- in post-smooth korakov
% maxLevel  … skupno št. nivojev

  % 1) PRED-Smoothing
  for i=1:nu1
    u = relaxGaussSeidel(A{level}, u, f, omega);
  end

  if level == maxLevel
    % 2a) NAJGRUBIJI NIVO: neposredno (GS) popravilo
    for i=1:nu2
      u = relaxGaussSeidel(A{level}, u, f, omega);
    end

  else
    % 2b) IZRAČUN RESIDUALA
    r = f - A{level}*u;

    % 3) Restrikcija na grobo mrežo
    rc = R{level} * r;

    % 4) Rekurzivni F-cikel na nivoju+1
    uc0 = zeros(size(rc));
    uc  = Fcycle(level+1, uc0, rc, A, R, P, omega, nu1, nu2, maxLevel);

    % 5) Prolongacija napake in korekcija
    u = u + P{level} * uc;

    % 6) POST-Smoothing
    for i=1:nu2
      u = relaxGaussSeidel(A{level}, u, f, omega);
    end
  end
end
%%---------------------------------------------------------------------------%%
