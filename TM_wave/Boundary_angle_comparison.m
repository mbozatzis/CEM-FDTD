addpath("TM_wave/PML/");

clc 
clear

cylinder_options = [0, 0, 1, 0, 0]; % Test boundary conditions with no cylinder
simulation_options = [10, 10, 100, 100, 12, 10*10^9];  
boundary_case = "full";  
PML_options = [16, 2, 10^(-6)]; 

% Calculate the Fields with different boundaries
boundary = "No-boundary"; 
[Ez_nb, Hx_nb, Hy_nb] = CylinderScattering_AngleTest(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

boundary = "Mur-first-order";
[Ez_m1, Hx_m1, Hy_m1] = CylinderScattering_AngleTest(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

boundary = "Mur-second-order";
[Ez_m2, Hx_m2, Hy_m2] = CylinderScattering_AngleTest(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

boundary = "PML";
[Ez_p, Hx_p, Hy_p] = CylinderScattering_AngleTest(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

% Calculate reference field
boundary = "No-boundary";
simulation_options = [25, 25, 250, 250, 12, 10*10^9];
[Ez_ref, Hx_ref, Hy_ref] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

% Define test points
tp_1 = [80, 100]; % ~90deg (grazing)
tp_2 = [100, 71]; % ~30deg
tp_3 = [100, 50]; % ~45deg
tp_4 = [100, 13]; % ~60deg

tp_1_ref = [155, 126];
tp_2_ref = [175, 96];
tp_3_ref = [175, 75];
tp_4_ref = [175, 38];

% Initialize error
er_nb = zeros(4, 1);
er_m1 = zeros(4, 1);
er_m2 = zeros(4, 1);
er_p = zeros(4, 1);

% 90 deg
Et_nb = reshape(Ez_nb(tp_1(2), tp_1(1), :), [1, size(Ez_nb, 3)]); 
Et_m1 = reshape(Ez_m1(tp_1(2), tp_1(1), :), [1, size(Ez_m1, 3)]);
Et_m2 = reshape(Ez_m2(tp_1(2), tp_1(1), :), [1, size(Ez_m2, 3)]);
Et_p = reshape(Ez_p(tp_1(2), tp_1(1), :), [1, size(Ez_p, 3)]);

Et_ref = reshape(Ez_ref(tp_1_ref(2), tp_1_ref(1), :), [1, size(Ez_ref, 3)]);

% Plot the errors of the boundary over time
t = linspace(0, 10, size(Ez_ref, 3));
ref = max(Et_ref);
figure(1);
sgtitle("Relative errors for 90° (Grazing)");
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

er_nb(4) = mean(nonzeros(abs(Et_ref - Et_nb)/ref));
er_m1(4) = mean(nonzeros(abs(Et_ref - Et_m1)/ref));
er_m2(4) = mean(nonzeros(abs(Et_ref - Et_m2)/ref));
er_p(4) = mean(nonzeros(abs(Et_ref - Et_p)/ref));


% 30 deg
Et_nb = reshape(Ez_nb(tp_2(2), tp_2(1), :), [1, size(Ez_nb, 3)]); 
Et_m1 = reshape(Ez_m1(tp_2(2), tp_2(1), :), [1, size(Ez_m1, 3)]);
Et_m2 = reshape(Ez_m2(tp_2(2), tp_2(1), :), [1, size(Ez_m2, 3)]);
Et_p = reshape(Ez_p(tp_2(2), tp_2(1), :), [1, size(Ez_p, 3)]);

Et_ref = reshape(Ez_ref(tp_2_ref(2), tp_2_ref(1), :), [1, size(Ez_ref, 3)]);

% Plot the errors of the boundary over time
t = linspace(0, 10, size(Ez_ref, 3));
ref = max(Et_ref);
figure(2);
sgtitle("Relative errors for 30°");
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

er_nb(1) = mean(nonzeros(abs(Et_ref - Et_nb)/ref));
er_m1(1) = mean(nonzeros(abs(Et_ref - Et_m1)/ref));
er_m2(1) = mean(nonzeros(abs(Et_ref - Et_m2)/ref));
er_p(1) = mean(nonzeros(abs(Et_ref - Et_p)/ref));


% 45 deg
Et_nb = reshape(Ez_nb(tp_3(2), tp_3(1), :), [1, size(Ez_nb, 3)]); 
Et_m1 = reshape(Ez_m1(tp_3(2), tp_3(1), :), [1, size(Ez_m1, 3)]);
Et_m2 = reshape(Ez_m2(tp_3(2), tp_3(1), :), [1, size(Ez_m2, 3)]);
Et_p = reshape(Ez_p(tp_3(2), tp_3(1), :), [1, size(Ez_p, 3)]);

Et_ref = reshape(Ez_ref(tp_3_ref(2), tp_3_ref(1), :), [1, size(Ez_ref, 3)]);

% Plot the errors of the boundary over time
t = linspace(0, 10, size(Ez_ref, 3));
ref = max(Et_ref);
figure(3);
sgtitle("Relative errors for 45°");
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

er_nb(2) = mean(nonzeros(abs(Et_ref - Et_nb)/ref));
er_m1(2) = mean(nonzeros(abs(Et_ref - Et_m1)/ref));
er_m2(2) = mean(nonzeros(abs(Et_ref - Et_m2)/ref));
er_p(2) = mean(nonzeros(abs(Et_ref - Et_p)/ref));

% 60 deg
Et_nb = reshape(Ez_nb(tp_4(2), tp_4(1), :), [1, size(Ez_nb, 3)]); 
Et_m1 = reshape(Ez_m1(tp_4(2), tp_4(1), :), [1, size(Ez_m1, 3)]);
Et_m2 = reshape(Ez_m2(tp_4(2), tp_4(1), :), [1, size(Ez_m2, 3)]);
Et_p = reshape(Ez_p(tp_4(2), tp_4(1), :), [1, size(Ez_p, 3)]);

Et_ref = reshape(Ez_ref(tp_4_ref(2), tp_4_ref(1), :), [1, size(Ez_ref, 3)]);

% Plot the errors of the boundary over time
t = linspace(0, 10, size(Ez_ref, 3));
ref = max(Et_ref);
figure(4);
sgtitle("Relative errors for 60°");
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

er_nb(3) = mean(nonzeros(abs(Et_ref - Et_nb)/ref));
er_m1(3) = mean(nonzeros(abs(Et_ref - Et_m1)/ref));
er_m2(3) = mean(nonzeros(abs(Et_ref - Et_m2)/ref));
er_p(3) = mean(nonzeros(abs(Et_ref - Et_p)/ref));


% Error versus Angle
figure(6);
angles = [30, 45, 60, 90];
plot(angles, er_nb);
hold on;
plot(angles, er_m1);
plot(angles, er_m2);
plot(angles, er_p);
legend('PEC', 'Mur 1st', 'Mur 2nd', 'PML');
title("Mean Relative Errors");
xlabel("Incident Angle (deg)");

