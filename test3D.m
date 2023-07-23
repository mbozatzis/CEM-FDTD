clc;
clear;

% Extract options
Xm_nl = 5;
Ym_nl = 5;
Zm_nl = 5;
N_x = 20;
N_y = 20;
N_z = 20;
Tn = 5;
f0 = 2*10^9;

% Useful constants
e0=8.85418781762e-12;
m0=1.256637061436e-06;
c0 = 3*10^8;
lambda0 = c0/f0;

% Simulation Characteristics
X_m = Xm_nl*lambda0;
Y_m = Ym_nl*lambda0;
Z_m = Zm_nl*lambda0;
dx = X_m/N_x;
dy = Y_m/N_y;
dz = Z_m/N_z;

Tmax = Tn/f0;
dt = 0.9*dx/(sqrt(3)*c0);

% Field matrices
Ex = zeros(N_x+1, N_y+1, N_z+1);
Ey = zeros(N_x+1, N_y+1, N_z+1);
Ez = zeros(N_x+1, N_y+1, N_z+1);

Hx = zeros(N_x+1, N_y, N_z);
Hy = zeros(N_x, N_y+1, N_z);
Hz = zeros(N_x, N_y, N_z+1);

% e = e0 * ones(N_x, N_y, N_z);

% Constants for field updating
Ce = dt/e0;
Cm = dt/m0;

[X, Y, Z] = meshgrid(1:size(Ex, 2), 1:size(Ex, 1), 1:size(Ex, 3));


for t = 0:dt:Tmax
    % Update Electric Field
    % Ex
    for i = 2:N_x
        for j = 2:N_y
            for mk = 2:N_z
                Ex(i, j, mk) = Ex(i, j, mk) + (Ce/dy)*(Hz(i, j, mk) - ...
                    Hz(i, j - 1, mk) - Hy(i, j, mk) + Hy(i, j, mk-1));
            end
        end
    end

    % Ey
    for i = 2:N_x
        for j = 2:N_y
            for mk = 2:N_z
                Ey(i, j, mk) = Ey(i, j, mk) + (Ce/dy)*(Hx(i, j, mk) - ...
                    Hx(i, j, mk-1) - Hz(i, j, mk) + Hz(i-1, j, mk));
            end
        end
    end

    % Ez
    for i = 2:N_x
        for j = 2:N_y
            for mk = 2:N_z
                Ez(i, j, mk) = Ez(i, j, mk) + (Ce/dy)*(Hy(i, j, mk) - ...
                    Hy(i-1, j, mk) - Hx(i, j, mk) + Hx(i, j-1, mk));
            end
        end
    end

    % Update source
    Ez(floor(N_x/2), floor(N_y/2), floor(N_z/2)) = sin(2*pi*f0*t);

    % Update Magnetic Fields
    % Hx
    for i = 2:N_x
        for j = 2:N_y
            for mk = 2:N_z
                Hx(i, j, mk) = Hx(i, j, mk) + (Cm/dy)*(Ey(i, j, mk+1) - ...
                    Ey(i, j, mk) - Ez(i, j+1, mk) + Ez(i, j, mk));
            end
        end
    end

    % Hy
    for i = 2:N_x
        for j = 2:N_y
            for mk = 2:N_z
                Hy(i, j, mk) = Hy(i, j, mk) + (Cm/dy)*(Ez(i+1, j, mk) - ...
                    Ez(i, j, mk) - Ex(i, j, mk+1) + Ex(i, j, mk));
            end
        end
    end

    % Hz
    for i = 2:N_x
        for j = 2:N_y
            for mk = 2:N_z
                Hz(i, j, mk) = Hz(i, j, mk) + (Cm/dy)*(Ex(i, j+1, mk) - ...
                    Ex(i, j, mk) - Ey(i+1, j, mk) + Ey(i, j, mk));
            end
        end
    end
end



