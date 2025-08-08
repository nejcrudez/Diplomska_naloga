function u = twoGridCycle(u, f, A1, A2, R1, P1, omega, nu1, nu2)
% En two-grid V-cikel: 
%  - nu1 GS pred-smooth na fine
%  - restrikcija reziduala
%  - solve na coarse (GS za nu1+nu2 korakov)
%  - prolongacija in korekcija
%  - nu2 GS post-smooth na fine

  % 1) pred-smooth
  for k = 1:nu1
    u = relaxGaussSeidel(A1, u, f, omega);
  end

  % 2) restrikcija reziduala
  r  = f - A1*u;
  rc = R1 * r;

  % 3) coarse solve (GS za nu1+nu2 korakov)
  uc = zeros(size(rc));
  for k = 1:(nu1+nu2)
    uc = relaxGaussSeidel(A2, uc, rc, omega);
  end

  % 4) prolongacija in popravilo
  u = u + P1 * uc;

  % 5) post-smooth
  for k = 1:nu2
    u = relaxGaussSeidel(A1, u, f, omega);
  end
end
