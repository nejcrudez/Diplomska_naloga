
n = 9; % velikost mreže (vedno moramo izbrati tako, da je n-1 sodo število)
A = poisson_stencil2D(n);  % matrika za Poissonovo enačbo
f = ones((n-1)^2, 1);  % desna stran enačbe (vse enice)
u = zeros(size(f));  % začetna rešitev
omega = 1;  % parameter Gauss-Seidela

z_min = 0;   % zgornja in spodnja meja za približke
z_max = 0.07;

% 100 korakov GS
for i = 1:100
    u = relaxGaussSeidel(A, u, f, omega);

    % Vizualizacija po vsakih 10 korakih

    if mod(i, 10) == 0
        U = reshape(u, n-1, n-1);  % ta del pretvori rešitev v 2D
        figure;
        surf(U);
        title(['Rešitev po ', num2str(i), ' korakih']);
        xlabel('x'); ylabel('y'); zlabel('u(x, y)');
        zlim([z_min, z_max]);     % višino fiksiramo
        clim([z_min, z_max]);
        colorbar;
        pause(0.1);
    end
end

% Izpis končne rešitve
disp('Končna rešitev (u):');
disp(u);

% Rezidual
r = f - A * u;

disp('Rezidual (r):');
disp(r);

% Prenos reziduala na grobo mrežo -
R = restrictionFW2D(n);  
r_coarse = R * r;
disp('Rezidual na grobi mreži (r_coarse):');
disp(r_coarse);
surf(reshape(r_coarse, 4, 4))
title('Rezidual na grobi mreži')
colorbar

T = interpolation2D((n-1)/2);  % interpolacijska matrika za grobo mrežo velikosti (n-1)/2
e_coarse = r_coarse;  % recimo, da je napaka približno enaka rezidualu (za poenostavitev problema)
e_fine = T * e_coarse;

% Popravek rešitve
u_corrected = u + e_fine;

% Vizualizacija popravka in popravljene rešitve
figure;
surf(reshape(e_fine, n-1, n-1));
title('Interpoliran popravek');
colorbar;

figure;
surf(reshape(u_corrected, n-1, n-1));
title('Popravljena rešitev po interpolaciji');
xlabel('x'); ylabel('y'); zlabel('u(x, y)');
colorbar;

 %Lahko tudi ponovno izračunamo rezidual po popravku (ni pa nujno):
r_corrected = f - A * u_corrected;
disp('Rezidual po popravku (r\_corrected):');
disp(r_corrected);