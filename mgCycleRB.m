function Unew = mgCycleRB(level, U, f, A, R, P, omega, nu1, nu2, gamma, maxLevel, h)
    % Recursive gamma-cycle with GS-RB on finest, GS on coarse
    % 1) pre-smoothing
    for i = 1:nu1
        if level==1
            U = relaxGaussSeidelRB(U, f{1}, omega, h{1});
        else
            uvec = relaxGaussSeidel(A{level}, U(:), f{level}(:), omega);
            U    = reshape(uvec, size(U));
        end
    end
    % 2) coarsest: only post-smoothing
    if level==maxLevel
        for i = 1:nu2
            if level==1
                U = relaxGaussSeidelRB(U, f{1}, omega, h{1});
            else
                uvec = relaxGaussSeidel(A{level}, U(:), f{level}(:), omega);
                U    = reshape(uvec, size(U));
            end
        end
        Unew = U; return;
    end
    % 3) compute residual & restrict
    if level==1
        D = f{1} - applyLaplace2D(U, h{1});
        r = D(:);
    else
        r = f{level}(:) - A{level}*U(:);
    end
    r = reshape(r, [], 1);
    rc = R{level} * r;
    
    
    % 4) recursion gamma-times
    ec = zeros(size(rc));
    for j = 1:gamma
        E  = reshape(ec, size(f{level+1}));  % Popravljeno
        ec = mgCycleRB(level+1, E, f, A, R, P, omega, nu1, nu2, gamma, maxLevel, h);
        ec = ec(:);
    end
    % 5) prolong & correct
    U = U + reshape(P{level}*ec, size(U));
    % 6) post-smoothing
    for i = 1:nu2
        if level==1
            U = relaxGaussSeidelRB(U, f{1}, omega, h{1});
        else
            uvec = relaxGaussSeidel(A{level}, U(:), f{level}(:), omega);
            U    = reshape(uvec, size(U));
        end
    end
    Unew = U;
end
