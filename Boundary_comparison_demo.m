% Before applying this code, the time-step dt was modified to be the same
% for every dx to have the same time reference for both cases.
% More specifically it was set to dt = 0.9*(lambda0/16)/(sqrt(2)*c0);
% in order to have stability in both cases (10x10, 16x16)

clc
clear

c0 = 3*10^8;
f0 = 10*10^9;
lambda0 = c0/f0;

Xm_nl = 10;
Ym_nl = 10;
dx = lambda0/Xm_nl;
n_lambda = lambda0/dx;

Xm_nl_ref = 16;
Ym_nl_ref = 16;
dx_ref = lambda0/Xm_nl_ref;
n_lambda_ref = lambda0/dx_ref;

cylinder_options = [0, 0, 1, 0, 0]; % Test boundary conditions with no cylinder
simulation_options = [10, 10, 10, 10*10^9];  
boundary_case = "full";  
PML_options = [8, 2, 10^(-6)]; 


boundary = "No-boundary"; 
[Ez_nb, Hx_nb, Hy_nb] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

boundary = "Mur-first-order";
[Ez_m1, Hx_m1, Hy_m1] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

boundary = "Mur-second-order";
[Ez_m2, Hx_m2, Hy_m2] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

boundary = "PML";
[Ez_p, Hx_p, Hy_p] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

figure(1);
subplot(2, 2, 1);
pcolor(Ez_nb(:, :, end)');
shading interp;
colormap default;
colorbar;
title('Zero boundary conditions');
xlabel('X-axis');
ylabel('Y-axis');

subplot(2, 2, 2);
pcolor(Ez_m1(:, :, end)');
shading interp;
colormap default;
colorbar;
title('Mur first order');
xlabel('X-axis');
ylabel('Y-axis');

subplot(2, 2, 3);
pcolor(Ez_m2(:, :, end)');
shading interp;
colormap default;
colorbar;
title('Mur second order');
xlabel('X-axis');
ylabel('Y-axis');

subplot(2, 2, 4);
pcolor(Ez_p(:, :, end)');
shading interp;
colormap default;
colorbar;
title('PML');
xlabel('X-axis');
ylabel('Y-axis');

% Reference
boundary = "No-boundary";
simulation_options = [16, 16, 10, 10*10^9];
[Ez_ref, Hx_ref, Hy_ref] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

test_point = [n_lambda*Xm_nl/2, 1]; % The middle of the bottom boundary
ref_test_point = [n_lambda_ref*Xm_nl_ref/2, ceil(n_lambda_ref*Ym_nl_ref/2 - (10/16)*n_lambda_ref*Ym_nl_ref/2)+1];

Et_nb = reshape(Ez_nb(test_point(1), test_point(2), :), [1, length(Ez_nb)]); 
Et_m1 = reshape(Ez_m1(test_point(1), test_point(2), :), [1, length(Ez_m1)]);
Et_m2 = reshape(Ez_m2(test_point(1), test_point(2), :), [1, length(Ez_m2)]);
Et_p = reshape(Ez_p(test_point(1), test_point(2), :), [1, length(Ez_p)]);

Et_ref = reshape(Ez_ref(ref_test_point(1), ref_test_point(2), :), [1, size(Ez_ref, 3)]);

ref = max(Et_ref);
figure(2);
sgtitle("Relative errors for a point in the middle of the bottom boundary");
subplot(2, 2, 1);
plot(abs(Et_ref - Et_nb)/ref);
title("No-boundary error");

subplot(2, 2, 2);
plot(abs(Et_ref - Et_m1)/ref);
title("Mur first order error");

subplot(2, 2, 3);
plot(abs(Et_ref - Et_m2)/ref);
title("Mur second order error");

subplot(2, 2, 4);
plot(abs(Et_ref - Et_p)/ref);
title("PML error");


