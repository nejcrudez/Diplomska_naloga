n = 9;
A = poisson_stencil2D(n);       % Matrika za Poissonovo enačbo
f = ones((n-1)^2, 1);           % Desna stran
u = zeros(size(f));             % Začetna rešitev
omega = 1.5;                      % parameter Gauss-Seidla (za omego, večje od 1, dobimo hitrejšo konvergenco) 

z_min = 0;
z_max = 0.07;

residuals = zeros(100, 1);

for i = 1:100
    u = relaxGaussSeidel(A, u, f, omega);

    % Izračunaj residual
    r = f - A * u;
    residuals(i) = norm(r, 2);

    % 2D prikaz rešitve vsakih 10 korakov
    if mod(i, 10) == 0
        U = reshape(u, n-1, n-1);
        figure;
        imagesc(U);
        colormap(parula);       % Izbira barvne lestvice
        colorbar;
        caxis([z_min z_max]);   % Fiksna barvna skala
        axis equal tight;
        title(['Stanje po ', num2str(i), ' korakih']);
        xlabel('x'); ylabel('y');
        pause(0.1);
    end
end

% Prikaz konvergence
figure;
semilogy(1:100, residuals, 'b-', 'LineWidth', 2);
xlabel('Iteracija');
ylabel('||f - A*u||_2');
title('Konvergenca Gauss-Seidla');
grid on;
