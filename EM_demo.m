
% To-dos
% CEM: 
% - Start making main (plots comparing boundaries, comparing different
% cylinders, time, PML options, put cylinder in different places, maybe an
% angle diagram for the boundary conditions, show PML plot, show stability and
% make it unstable in purpose, velocity plot (kdx and angle), maybe frequency?)
% - report with theoretical calculations
%
% Scattering:
% - TF/SF Formulation implementation
% - PML at the corners
% - Analytical solution implementation
% - Simulation with different cylinders at different positions and angles
% - Add complex materials

clc
clear


% Part 1 - No boundaries
cylinder_options = [3, 1, 3.4, 1.2]; % [x0, R, e_r, sigma_r]
simulation_options = [10, 10, 12, 10*10^9]; % [X0, Y0, Tn, f0]
boundary = "No-boundary"; % Type of boundary condition
boundary_case = "full"; % Not in use in case of no-boundaries
PML_options = [8, 2, 10^(-6)]; % Only in use when boundary is PML

[Ez_nb, Hx_nb, Hy_nb] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

Tn = 12; % Total simulation time: 12 T0
i_max = size(Ez_nb, 3);
ind_t1 = ceil((3/12)*i_max); % t1 = 3T0
ind_t2 = ceil((10/12)*i_max); % t2 = 10T0
ind_t3 = i_max; % t3 = 12T0 = TnT0

figure(1);
sgtitle('Electric Field with zero boundary conditions');

subplot(2, 2, 1);
pcolor(Ez_nb(:, :, ind_t1));
shading interp;
colormap default;
colorbar;
title('Ez field at time: t1 = 3T0');
xlabel('X-axis');
ylabel('Y-axis');

subplot(2, 2, 2);
pcolor(Ez_nb(:, :, ind_t2));
shading interp;
colormap default;
colorbar;
title('Ez field at time: t2 = 10T0');
xlabel('X-axis');
ylabel('Y-axis');

subplot(2, 2, 3);
pcolor(Ez_nb(:, :, ind_t3));
shading interp;
colormap default;
colorbar;
title('Ez field at time: t3 = 12T0');
xlabel('X-axis');
ylabel('Y-axis');


% Part 2 - Mur boundaries
cylinder_options = [3, 1, 3.4, 1.2]; % [x0, R, e_r, sigma_r]
simulation_options = [10, 10, 12, 10*10^9]; % [X0, Y0, Tn, f0]
boundary = "Mur-first-order"; % Type of boundary condition
boundary_case = "not-full"; % boundaries only in the left side
PML_options = [8, 2, 10^(-6)]; % Only in use when boundary is PML

[Ez_m1, Hx_m1, Hy_m1] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

Tn = 12; % Total simulation time: 12 T0
i_max = size(Ez_m1, 3);
ind_t1 = ceil((6/12)*i_max); % t1 = 6T0
ind_t2 = i_max; % t2 = 12T0 = TnT0

figure(2);
sgtitle('Electric Field with first order Mur boundary condition on left side');

subplot(1, 2, 1);
pcolor(Ez_m1(:, :, ind_t1));
shading interp;
colormap default;
colorbar;
title('Ez field at time: t1 = 6T0');
xlabel('X-axis');
ylabel('Y-axis');

subplot(1, 2, 2);
pcolor(Ez_m1(:, :, ind_t2));
shading interp;
colormap default;
colorbar;
title('Ez field at time: t2 = 12T0');
xlabel('X-axis');
ylabel('Y-axis');

boundary = "Mur-second-order";
[Ez_m2, Hx_m2, Hy_m2] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

figure(3);
sgtitle('Electric Field with second order Mur boundary condition on left side');

subplot(1, 2, 1);
pcolor(Ez_m2(:, :, ind_t1));
shading interp;
colormap default;
colorbar;
title('Ez field at time: t1 = 6T0');
xlabel('X-axis');
ylabel('Y-axis');

subplot(1, 2, 2);
pcolor(Ez_m2(:, :, ind_t2));
shading interp;
colormap default;
colorbar;
title('Ez field at time: t2 = 12T0');
xlabel('X-axis');
ylabel('Y-axis');


% Part 3 - PML
cylinder_options = [3, 1, 3.4, 1.2]; % [x0, R, e_r, sigma_r]
simulation_options = [10, 10, 12, 10*10^9]; % [X0, Y0, Tn, f0]
boundary = "PML"; % Type of boundary condition
boundary_case = "not-full"; % boundaries only in the left side
PML_options = [8, 2, 10^(-6)]; % [Npml, power, R]

[Ez_p, Hx_p, Hy_p] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

Tn = 12; % Total simulation time: 12 T0
i_max = size(Ez_p, 3);
ind_t1 = ceil((6/12)*i_max); % t1 = 6T0
ind_t2 = i_max; % t2 = 12T0 = TnT0

figure(4);
sgtitle('Electric Field with PML on left side');

subplot(1, 2, 1);
pcolor(Ez_p(:, :, ind_t1));
shading interp;
colormap default;
colorbar;
title('Ez field at time: t1 = 6T0');
xlabel('X-axis');
ylabel('Y-axis');

subplot(1, 2, 2);
pcolor(Ez_p(:, :, ind_t2));
shading interp;
colormap default;
colorbar;
title('Ez field at time: t2 = 12T0');
xlabel('X-axis');
ylabel('Y-axis');

% Part 4 - Electric field at points of grid over time
c0 = 3*10^8;
f0 = 10*10^9;
lambda0 = c0/f0;
Xm_nl = 10;
Ym_nl = 10;
dx = lambda0/Xm_nl;
n_lambda = lambda0/dx;

x1 = n_lambda;
x2 = n_lambda;
y1 = ceil(Ym_nl/2);
y2 = ceil(Ym_nl/2) + n_lambda;

figure(5);
sgtitle('Electric Field at points of grid over time');

subplot(4, 2, 1);
plot(reshape(Ez_nb(x1, y1, :), [1, length(Ez_nb)]));
title('Point x1 with no boundary conditions');
ylabel('Ez field');
xlabel('time');

subplot(4, 2, 2);
plot(reshape(Ez_nb(x2, y2, :), [1, length(Ez_nb)]));
title('Point x2 with no boundary conditions');
ylabel('Ez field');
xlabel('time');

subplot(4, 2, 3);
plot(reshape(Ez_m1(x1, y1, :), [1, length(Ez_m1)]));
title('Point x1 with first order Mur');
ylabel('Ez field');
xlabel('time');

subplot(4, 2, 4);
plot(reshape(Ez_m1(x2, y2, :), [1, length(Ez_m1)]));
title('Point x2 with first order Mur');
ylabel('Ez field');
xlabel('time');

subplot(4, 2, 5);
plot(reshape(Ez_m2(x1, y1, :), [1, length(Ez_m2)]));
title('Point x1 with second order Mur');
ylabel('Ez field');
xlabel('time');

subplot(4, 2, 6);
plot(reshape(Ez_m2(x2, y2, :), [1, length(Ez_m2)]));
title('Point x2 with second order Mur');
ylabel('Ez field');
xlabel('time');

subplot(4, 2, 7);
plot(reshape(Ez_p(x1, y1, :), [1, length(Ez_p)]));
title('Point x1 PML');
ylabel('Ez field');
xlabel('time');

subplot(4, 2, 8);
plot(reshape(Ez_p(x2, y2, :), [1, length(Ez_p)]));
title('Point x2 with PML');
ylabel('Ez field');
xlabel('time');




% writerObj = VideoWriter('ElectricField.avi');
% writerObj.FrameRate = 25; 
% open(writerObj);
% 
% for k = 1:size(Ez,3)
%     pcolor(Ez(:, :, k));
%     shading interp;
%     colormap gray;
%     
%     title(['Contour plot at slice ' num2str(k)]);
%     xlabel('X-axis');
%     ylabel('Y-axis');
%     
%     % Capture the frame and add it to the video
%     currFrame = getframe(gcf);
%     writeVideo(writerObj, currFrame);
% end
% 
% close(writerObj);






