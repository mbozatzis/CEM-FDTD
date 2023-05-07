clc
clear

c0 = 3*10^8;
f0 = 10*10^9;
lambda0 = c0/f0;
Xm_nl = 10;
Ym_nl = 10;
dx = lambda0/Xm_nl;
n_lambda = lambda0/dx;

theta = linspace(0, 2*pi);
x_c = (Xm_nl/2-3)*n_lambda + n_lambda*cos(theta);
y_c = (Ym_nl/2)*n_lambda + n_lambda*sin(theta);


simulation_options = [10, 10, 12, 10*10^9]; % [X0, Y0, Tn, f0]
boundary = "Mur-second-order"; % Type of boundary condition
boundary_case = "full"; % Second order Mur at all boundaries 
PML_options = [8, 2, 10^(-6)]; % Only in use when boundary is PML

% Nearly perfect conductor
cylinder_options = [3, 1, 3.4, 100, 0];
[Ez, Hx, Hy] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

i_max = size(Ez, 3);

figure(1);
pcolor(Ez(:, :, i_max)');
hold on;
plot(x_c, y_c, 'r', 'LineWidth', 1);
title('Perfectly conducting scatterer');
shading interp;
colormap default;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');

% Dielectric cylinder
cylinder_options = [3, 1, 3.4, 0, 0];
[Ez, Hx, Hy] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

i_max = size(Ez, 3);

figure(2);
pcolor(Ez(:, :, i_max)');
hold on;
plot(x_c, y_c, 'r', 'LineWidth', 1);
title('Dielectric scatterer');
shading interp;
colormap default;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');

% Scatterer at different place
cylinder_options = [3, 1, 3.4, 1.2, 3];
[Ez, Hx, Hy] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

x_c = (Xm_nl/2-3)*n_lambda + n_lambda*cos(theta);
y_c = (Ym_nl/2-3)*n_lambda + n_lambda*sin(theta);

i_max = size(Ez, 3);

figure(3);
sgtitle('Cylinder at different place');

subplot(2, 2, 1);
pcolor(Ez(:, :, i_max)');
hold on;
plot(x_c, y_c, 'r', 'LineWidth', 1);
title('Default scatterer');
shading interp;
colormap default;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');

cylinder_options = [3, 1, 3.4, 100, 3];
[Ez, Hx, Hy] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

subplot(2, 2, 2);
pcolor(Ez(:, :, i_max)');
hold on;
plot(x_c, y_c, 'r', 'LineWidth', 1);
title('Perfectly conducting scatterer');
shading interp;
colormap default;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');

cylinder_options = [3, 1, 3.4, 0, 3];
[Ez, Hx, Hy] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

subplot(2, 2, 3);
pcolor(Ez(:, :, i_max)');
hold on;
plot(x_c, y_c, 'r', 'LineWidth', 1);
title('Dielectric scatterer');
shading interp;
colormap default;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');

subplot(2, 2, 4);
pcolor(Hy(:, :, i_max)');
hold on;
plot(x_c, y_c, 'r', 'LineWidth', 1);
title('Dielectric scatterer - Magnetic Field Hy');
shading interp;
colormap default;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');