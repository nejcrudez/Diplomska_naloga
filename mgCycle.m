function u = mgCycle(level, u, f, A, R, P, omega, nu1, nu2, gamma, maxLevel)
% MGCYCLE  Generaliziran multigrid cikel s parametrom gamma.
% level     ... trenutni nivo (1=finest)
% u, f      ... vektor rešitve in desna stran na tem nivoju
% A, R, P   ... celice z operatorji za vse nivoje
% omega     ... SOR parameter
% nu1,nu2   ... št. pred- in post-smooth GS korakov
% gamma     ... rekurzivna globina (1=V, 2=W, ...)
% maxLevel  ... število nivojev

  % 1) pred-glajenje
  for i = 1:nu1
    u = relaxGaussSeidel(A{level}, u, f, omega);
  end

  if level == maxLevel
    % coarsest nivo: samo post-smoothing
    for i = 1:nu2
      u = relaxGaussSeidel(A{level}, u, f, omega);
    end
  else
    % 2) izračun reziduala
    r = f - A{level}*u;
    % 3) restrikcija
    rc = R{level}*r;
    % 4) rekurzija gamma-krat
    ec = zeros(size(rc));
    for j = 1:gamma
      ec = mgCycle(level+1, ec, rc, A, R, P, omega, nu1, nu2, gamma, maxLevel);
    end
    % 5) prolongacija + korekcija
    u = u + P{level} * ec;
    % 6) post-glajenje
    for i = 1:nu2
      u = relaxGaussSeidel(A{level}, u, f, omega);
    end
  end
end
