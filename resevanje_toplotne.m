%% 1) Začetne nastavitve
clear; clc;

n        = 32;
h        = 1/n;
N        = (n-1)^2;
T        = 0.1;

tau_list = [1/1000, 1/5000, 1/20000];   % primer treh tau vrednosti

%% 2) Operator Laplace
% Predpostavka: poisson_stencil2D(m) vrne sparse (n-1)^2 x (n-1)^2 matriko
% z Dirichlet robnimi pogoji na enotskem kvadratu, skladno s tvojim solverjem.
Lop = @(m) poisson_stencil2D(m);

%% 3) Začetni pogoj (disk v središču)
[x,y] = meshgrid(h:h:1-h, h:h:1-h);
U0 = double((x - 0.5).^2 + (y - 0.5).^2 <= 0.1);  % binarna "kapljica"
figure('Color','w');
surf(x, y, U0, 'EdgeColor','none');
caxis([0,1]); colormap(jet); view(2); colorbar;
title('t = 0.000');
drawnow;

%% 3b) Evolucija z “direct solverjem” (Backward Euler za izbrani tau)
tau = 1/5000;                         % izberi en tau za prikaz
A   = -speye((n-1)^2) - tau*Lop(n);   % A = -(I + tau*L), skladno s tvojim zapisom
u_direct = U0(:);                     % pravilni začetni vektor

numSteps_direct = round(T/tau);
for t = 1:numSteps_direct
    u_direct = A \ u_direct;          % reši A * u^{k+1} = u^{k}
    if mod(t, 50) == 0 || t == numSteps_direct
        Uplot = reshape(u_direct, n-1, n-1);
        figure('Color','w');
        surf(x, y, Uplot, 'EdgeColor','none');
        caxis([0,1]); colormap(jet);
        view(2); colorbar; 
        title(sprintf('Direct solver, t = %.4f (tau = %.5f)', t*tau, tau));
        drawnow;
    end
end

%% 4) Zanka čez različne tau za večmrežno metodo (MG)
tol     = 1e-8;
maxIter = 10;
results = zeros(length(tau_list), 3); % stolpci: [avg, min, max] št. MG iteracij na časovni korak

for k = 1:length(tau_list)
    tau = tau_list(k);
    numSteps = round(T/tau);

    % A(m) naj bo konsistenten z zgornjim direktnim solverjem
    Aop = @(m) -speye((m-1)^2) - tau*Lop(m);

    % začetni vektor
    u_mg = U0(:);
    MG_iters_per_step = zeros(numSteps, 1);

    for t = 1:numSteps
        f = u_mg;   % desna stran: A * u^{k+1} = u^{k}

        % MG solver iteracije do tolerance ali maxIter
        iter = 0;
        res  = inf;

        while (res > tol) && (iter < maxIter)
            % Predpostavka: multigrid_cycle vrne posodobljen približek rešitve
            % Signatura: multigrid_cycle(n, level, u, Aop, f, nu1, nu2, 'SmoothingMethod','GaussSeidel')
            u_mg = multigrid_cycle(n, 1, u_mg, Aop, f, 3, 3, 'SmoothingMethod','GaussSeidel');

            % ostanek (rezidual)
            res = norm(f - Aop(n) * u_mg);
            iter = iter + 1;
        end

        MG_iters_per_step(t) = iter;

        if mod(t, 50) == 0 || t == numSteps
            Uplot = reshape(u_mg, n-1, n-1);
            figure('Color','w');
            surf(x, y, Uplot, 'EdgeColor','none');
            caxis([0,1]); colormap(jet);
            view(2); colorbar; 
            % --- poenostavljen naslov brez dodatkov ---
            title(sprintf('Večmrežna metoda, tau = %.5f, t = %.4f', tau, t*tau));
            drawnow;
            fprintf('tau = %.5f, t = %.4f: MG iter = %d, res = %.3e\n', tau, t*tau, iter, res);
        end
    end

    % shranimo statistiko za ta tau
    results(k, :) = [mean(MG_iters_per_step), min(MG_iters_per_step), max(MG_iters_per_step)];
end

%% 5) Izpis tabele s statistiko MG
Tau        = tau_list(:);
Povprecje  = results(:,1);
Min        = results(:,2);
Max        = results(:,3);

T_table = table(Tau, Povprecje, Min, Max);
disp('Povzetek večmrežnih iteracij na časovni korak za različne tau:');
disp(T_table);
