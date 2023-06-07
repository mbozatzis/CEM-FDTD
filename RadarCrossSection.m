addpath("TFSF/");

clc
clear

cylinder_options_c = [0, 1, 1, 100, 0]; % [x0, R, e_r, sigma_r, y0]
cylinder_options_d = [0, 1, 1.6, 0, 0]; % [x0, R, e_r, sigma_r, y0]
simulation_options = [10, 10, 100, 100, 10, 10*10^9]; % [X0, Y0, N_x, N_y Tn, f0]
PML_options = [16, 2, 10^(-6)]; % [Npml, power, R]
TFSF_options = [40, pi/2]; % [Ntfsf, fi]

% With conducting cylinder
[eez_c, hhx_c, hhy_c] = TFSF_formulation(cylinder_options_c, simulation_options, ...
    PML_options, TFSF_options);

% With dielectric cylinder
[eez_d, hhx_d, hhy_d] = TFSF_formulation(cylinder_options_d, simulation_options, ...
    PML_options, TFSF_options);

c0 = 3*10^8;
f0 = 10*10^9;
T0 = 1/f0;
Tmax = 10*T0;
lambda0 = c0/f0;

Xm_nl = 10;
Ym_nl = 10;
N_x = 100;
N_y = 100;
X_m = Xm_nl*lambda0;
Y_m = Ym_nl*lambda0;
dx = X_m/N_x;
dy = Y_m/N_y;
n_lambda = lambda0/dx;
D = 2*n_lambda/10;
k0 = 2*pi/lambda0;

Ez_c = eez_c(:, :, end);
Ez_d = eez_d(:, :, end);

s_c = zeros();
s_d = zeros();
Theta = zeros();

% Define the range for x and y
xRange = linspace(-5, 5, 100);
yRange = linspace(-5, 5, 100);

% Create a grid of x and y values
[X, Y] = meshgrid(xRange, yRange);

% Calculate r and theta for each grid point
r = sqrt(X.^2 + Y.^2);
theta = atan2(Y, X);

l = 1;
for i = 1:100
    for j = 1:100
        s_c(l) = 2*pi*(r(i,j)*abs(Ez_c(i, j))^2);
        s_d(l) = 2*pi*(r(i,j)*abs(Ez_d(i, j))^2);
        Theta(l) = theta(i,j)*180/pi;
        l = l + 1;
    end
end

% Plot Bistable Radar crossection for Perfect Conductor
figure(1);
plot(Theta, s_c, '.');
title("Bistable Radar crossection");
xlabel("Angle (deg)");

% Bistable Radar crossection for Dielectric cylinder
figure(2);
plot(Theta, s_d, '.');
title("Bistable Radar crossection");
xlabel("Angle (deg)");


