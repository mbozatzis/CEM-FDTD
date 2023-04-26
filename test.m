clc
clear


e0=8.85418781762e-12;
m0=1.256637061436e-06;
c0 = 3*10^8;

f0 = 10*10^9;
lambda0 = c0/f0;

X_m = 10*lambda0;
Y_m = 10*lambda0;

% Cylinder position and characteristics
x0 = 3*lambda0;
R = lambda0;
e_r = 3.4;
sigma_r = 1.2;

dx = lambda0/10;
dy = lambda0/10;
N_x = ceil(X_m/dx);
N_y = ceil(Y_m/dy);

Tmax = 10/f0;
dt = 0.9*dx/(sqrt(2)*c0);
N_T = ceil(Tmax/dt);

% TM: Ez, Hx, Hy
Ez = zeros(N_x, N_y);
Hx = zeros(N_x, N_y);
Hy = zeros(N_x, N_y);

e = e0*zeros(N_x, N_y);
sigma = zeros(N_x, N_y);

% Create the cylinder
n_lambda = lambda0/dx;
i0 = N_x/2 - 3*n_lambda;
j0 = N_y/2;
r = n_lambda;
for i = 1:N_x
    for j = 1:N_y
        if sqrt((i-i0)^2 + (j-j0)^2) <= r
            e(i, j) = e_r*e0;
            sigma(i, j) = sigma_r;
        end
    end
end

% surf(e) % Visualize the space
% surf(sigma) % Visualize the space

