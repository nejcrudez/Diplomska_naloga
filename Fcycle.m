function u = Fcycle(k, u, f, A, R, P, omega, nu1, nu2, levels)
% En F-cikel multigrid (rekurzivno)

    % pred-smoothing
    for i = 1:nu1
        u = relaxGaussSeidel(A{k}, u, f, omega);
    end

    % rezidual in restrikcija
    r = f - A{k} * u;
    if k == levels
        % solve na najgrobnejši mreži (npr. nekaj GS korakov)
        e = zeros(size(r));
        for i = 1:(nu1 + nu2)
            e = relaxGaussSeidel(A{k}, e, r, omega);
        end
    else
        rc = R{k} * r;
        ec = zeros(size(rc));

        % REKURZIVNO: najprej navzdol (F-cikel zahteva 2 rekurziji)
        ec = Fcycle(k+1, ec, rc, A, R, P, omega, nu1, nu2, levels);
        ec = Fcycle(k+1, ec, rc, A, R, P, omega, nu1, nu2, levels);
        
        % korekcija
        e = P{k} * ec;
    end

    u = u + e;

    % post-smoothing
    for i = 1:nu2
        u = relaxGaussSeidel(A{k}, u, f, omega);
    end
end
