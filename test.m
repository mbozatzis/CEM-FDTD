clc
clear

cylinder_options = [0, 1, 3.4, 1.2, 0]; % [x0, R, e_r, sigma_r, y0]
simulation_options = [10, 10, 100, 100, 10, 10*10^9]; % [X0, Y0, N_x, N_y Tn, f0]
PML_options = [16, 2, 10^(-6)]; % [Npml, power, R]
TFSF_options = [8, pi/2]; % [Ntfsf, fi]

[eez, hhx, hhy] = TFSF_formulation(cylinder_options, simulation_options, ...
    PML_options, TFSF_options);
