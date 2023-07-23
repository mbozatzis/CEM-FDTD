addpath("PML/")
addpath("TFSF/");

clc
clear

cylinder_options = [0, 1, 1, 0, 0]; % [x0, R, e_r, sigma_r, y0]
cylinder_options_sc = [0, 1, 3.4, 1.2, 0];
simulation_options = [10, 10, 100, 100, 8, 10*10^9]; % [X0, Y0, N_x, N_y Tn, f0]
PML_options = [16, 2, 10^(-6)]; % [Npml, power, R]
TFSF_options_0 = [8, 0]; % [Ntfsf, fi]
TFSF_options_45 = [8, pi/4];
TFSF_options_90 = [8, pi/2];
TFSF_options_15 = [8, pi/12];

% Without scatterer
[eez_0, ~, ~] = TFSF_formulation(cylinder_options, simulation_options, ...
    PML_options, TFSF_options_0);
[eez_15, ~, ~] = TFSF_formulation(cylinder_options, simulation_options, ...
    PML_options, TFSF_options_15);
[eez_45, ~, ~] = TFSF_formulation(cylinder_options, simulation_options, ...
    PML_options, TFSF_options_45);
[eez_90, ~, ~] = TFSF_formulation(cylinder_options, simulation_options, ...
    PML_options, TFSF_options_90);

% With scatterer
[eez_0_sc, ~, ~] = TFSF_formulation(cylinder_options_sc, simulation_options, ...
    PML_options, TFSF_options_0);
[eez_15_sc, ~, ~] = TFSF_formulation(cylinder_options_sc, simulation_options, ...
    PML_options, TFSF_options_15);
[eez_45_sc, ~, ~] = TFSF_formulation(cylinder_options_sc, simulation_options, ...
    PML_options, TFSF_options_45);
[eez_90_sc, ~, ~] = TFSF_formulation(cylinder_options_sc, simulation_options, ...
    PML_options, TFSF_options_90);

% Time evolution of plane wave
i_max = size(eez_0, 3);
ind_t1 = ceil((3/8)*i_max); % t1 = 3T0
ind_t2 = ceil((5/8)*i_max); % t2 = 5T0

figure(1);
subplot(2, 1, 1);
pcolor(eez_0(:, :, ind_t1));
shading interp;
title("Plane Wave t_1 = 3T_0")

subplot(2, 1, 2);
pcolor(eez_0(:, :, ind_t2));
shading interp;
title("Plane Wave t_2 = 5T_0")


% Demo of plane wave generation
figure(2);
subplot(2, 2, 1);
pcolor(eez_0(:, :, end));
shading interp;
title("Plane Wave 0°")

subplot(2, 2, 2);
pcolor(eez_15(:, :, end));
shading interp;
title("Plane Wave 15°")

subplot(2, 2, 3);
pcolor(eez_45(:, :, end));
shading interp;
title("Plane Wave 45°")

subplot(2, 2, 4);
pcolor(eez_90(:, :, end));
shading interp;
title("Plane Wave 90°")

% Field at center over time
tp = [50, 50];
t = linspace(0, 8, size(eez_0, 3));

figure(2);
subplot(2, 2, 1);
plot(t, reshape(eez_0(tp(1), tp(2), :), [1, size(eez_0, 3)]));
title("Plane Wave 0°")

subplot(2, 2, 2);
plot(t, reshape(eez_15(tp(1), tp(2), :), [1, size(eez_15, 3)]));
title("Plane Wave 15°")

subplot(2, 2, 3);
plot(t, reshape(eez_45(tp(1), tp(2), :), [1, size(eez_45, 3)]));
title("Plane Wave 45°")

subplot(2, 2, 4);
plot(t, reshape(eez_90(tp(1), tp(2), :), [1, size(eez_90, 3)]));
title("Plane Wave 90°")


% Time evolution of plane wave to cylinder
i_max = size(eez_0_sc, 3);
ind_t1 = ceil((3/8)*i_max); % t1 = 3T0
ind_t2 = ceil((5/8)*i_max); % t2 = 5T0

figure(4);
subplot(2, 1, 1);
pcolor(eez_0_sc(:, :, ind_t1));
shading interp;
title("Plane Wave t_1 = 3T_0")

subplot(2, 1, 2);
pcolor(eez_0_sc(:, :, ind_t2));
shading interp;
title("Plane Wave t_2 = 5T_0")


% Demo of plane waves scattering to cylinder
figure(5);
subplot(2, 2, 1);
pcolor(eez_0_sc(:, :, end));
shading interp;
title("Plane Wave 0°")

subplot(2, 2, 2);
pcolor(eez_15_sc(:, :, end));
shading interp;
title("Plane Wave 15°")

subplot(2, 2, 3);
pcolor(eez_45_sc(:, :, end));
shading interp;
title("Plane Wave 45°")

subplot(2, 2, 4);
pcolor(eez_90_sc(:, :, end));
shading interp;
title("Plane Wave 90°")



