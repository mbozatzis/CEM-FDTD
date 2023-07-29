addpath("TM_wave/PML/");

clc
clear

cylinder_options = [0, 0, 1, 0, 0]; % Test boundary conditions with no cylinder
simulation_options = [10, 10, 100, 100, 10, 10*10^9];  
boundary_case = "full";  
PML_options = [16, 2, 10^(-6)]; 

% Calculate the Fields with different boundaries
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

% Plot the Ez Field with different boundaries
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

% Reference Field simulation
boundary = "No-boundary";
simulation_options = [16, 16, 160, 160, 10, 10*10^9];
[Ez_ref, Hx_ref, Hy_ref] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

% Add testpoint at the middle of the bottom boundary
test_point = [floor(size(Ez_nb, 1)/2), 1]; 

% Add the tespoint at the relevant position in the reference grid
ref_test_point = [floor(size(Ez_ref, 1)/2), floor(size(Ez_ref, 2)/2) - floor(size(Ez_nb, 2)/2) + 1];

% Save the field behaviour at the testpoint over time
Et_nb = reshape(Ez_nb(test_point(1), test_point(2), :), [1, length(Ez_nb)]); 
Et_m1 = reshape(Ez_m1(test_point(1), test_point(2), :), [1, length(Ez_m1)]);
Et_m2 = reshape(Ez_m2(test_point(1), test_point(2), :), [1, length(Ez_m2)]);
Et_p = reshape(Ez_p(test_point(1), test_point(2), :), [1, length(Ez_p)]);

Et_ref = reshape(Ez_ref(ref_test_point(1), ref_test_point(2), :), [1, size(Ez_ref, 3)]);

% Plot the errors of the boundary over time
t = linspace(0, 10, size(Ez_ref, 3));
ref = max(Et_ref);
figure(2);
sgtitle("Relative errors for different boundaries");
subplot(2, 2, 1);
plot(t, abs(Et_ref - Et_nb)/ref);
xlabel("Time (T_0)")
ylabel("(E_z - E_{ref})/E_{ref}")
title("No-boundary error");

subplot(2, 2, 2);
plot(t, abs(Et_ref - Et_m1)/ref);
xlabel("Time (T_0)")
ylabel("(E_z - E_{ref})/E_{ref}")
title("Mur first order error");

subplot(2, 2, 3);
plot(t, abs(Et_ref - Et_m2)/ref);
xlabel("Time (T_0)")
ylabel("(E_z - E_{ref})/E_{ref}")
title("Mur second order error");

subplot(2, 2, 4);
plot(t, abs(Et_ref - Et_p)/ref);
xlabel("Time (T_0)")
ylabel("(E_z - E_{ref})/E_{ref}")
title("PML error");




