
clear; clc;

%% 1) Začetne nastavitve
n        = 33;                 % število točk(notranjih je 32)
h        = 1/n;                % prostorski korak
A        = poisson_stencil2D(n);
tau      = 0.001;              % časovni korak
T        = 0.05;               % končni čas
omega    = 1;                  % Gauss–Seidel parameter
numSteps = round(T/tau);

%% 2) Definiramo toplo pego v središču
[x,y] = meshgrid(h:h:1-h, h:h:1-h);
U     = exp(-50*((x-0.5).^2 + (y-0.5).^2));
u     = U(:);

%% 3) Priprava figure
figure('Color','w');
surf(x, y, U, 'EdgeColor', 'none');
colormap(jet);
caxis([0, 0.07]);
view(2); colorbar;
title('t = 0.000');
drawnow;

%% 4) Časovna zanka
for t = 1:numSteps
    % napoved z eksplicitno Eulerjevo metodo
    f = u + tau*(A*u);
    
    % 4b) predglajenje: 10 korakov
    for k = 1:10
        u = relaxGaussSeidel(A, u, f, omega);
    end
    
    % 4c) rezidual
    r = f - A*u;
    
    % 4d) en V-cikel 
    R      = restrictionFW2D(n);
    Pmat   = interpolation2D((n-1)/2);
    rc     = R * r;
    ec     = rc;               % coarse solve ≈ residual
    e_fine = Pmat * ec;        % prolongacija
    u      = u + e_fine;       % popravek
    
    % 4e) vizualizacija na vsakem 5. koraku
    if mod(t, 5) == 0
        Uplot = reshape(u, n-1, n-1);
        figure(t)
        surf(x, y, Uplot, 'EdgeColor', 'none');
        title(sprintf('t = %.3f', t*tau));
        drawnow;
    end
end

%% 5) Končna rešitev
disp('Končna rešitev (u):');
Ufinal = reshape(u, n-1, n-1);
disp(Ufinal);

