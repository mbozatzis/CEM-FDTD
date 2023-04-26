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
Ez = zeros(N_x+1, N_y+1);
Hx = zeros(N_x+1, N_y);
Hy = zeros(N_x, N_y+1);

e = e0*ones(N_x+1, N_y+1);
sigma = zeros(N_x+1, N_y+1);

% Create the cylinder
n_lambda = lambda0/dx;
i0 = N_x/2 - 3*n_lambda;
j0 = N_y/2;
r = n_lambda;
for i = 1:N_x+1
    for j = 1:N_y+1
        if sqrt((i-i0)^2 + (j-j0)^2) <= r
            e(i, j) = e_r*e0;
            sigma(i, j) = sigma_r;
        end
    end
end

% surf(e) % Visualize the space
% surf(sigma) % Visualize the space

Ca = (2*e-dt*sigma)./(2*e+dt*sigma);
Cb = (dt./(e + dt * sigma/2))*(1/dx);
Da = dt/(m0*dy);
Db = dt/(m0*dx);

ee = zeros(N_x+1, N_y+1, N_T);
k = 1;
for t = 0:dt:Tmax
    for i = 2:N_x
        for j = 2:N_y
            Ez(i, j) = Ca(i, j)*Ez(i, j) + Cb(i, j)*(Hy(i, j)-Hy(i-1, j) ...
                + Hx(i,j-1) - Hx(i, j));
        end
    end

    Ez(N_x/2, N_y/2) = sin(2*pi*f0*t);

    for i = 2:N_x
        for j = 1:N_y
            Hx(i, j) = Hx(i, j) - Da*(Ez(i, j+1) - Ez(i, j));
        end
    end

    for i = 1:N_x
        for j = 2:N_y
            Hy(i, j) = Hy(i, j) + Db*(Ez(i+1, j) - Ez(i, j));
        end
    end
    
    surf(Ez);
    colormap default;
    drawnow;

    ee(:, :, k) = Ez;
    k = k +1;
end

