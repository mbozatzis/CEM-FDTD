clc
clear

% Define useful constants
f0 = 10*10^9;
e0=8.85418781762e-12;
m0=1.256637061436e-06;
c0 = 3*10^8;
lambda0 = c0/f0;

% Define simulation options
N_x = 100;
N_y = 100;
Xm_nl = 10;
Ym_nl = 10;
X_m = Xm_nl*lambda0;
Y_m = Ym_nl*lambda0;
dx = X_m/N_x;
dy = Y_m/N_y;
n_lambda = lambda0/dx;

i_center = floor(N_x/2);
j_center = floor(N_y/2);

n_num = 100;
Es = zeros(N_x, N_y);
Et = zeros(N_x, N_y);
Einc = zeros(N_x, N_y);
R = 1;
k0 = 2*pi;


% Define the range for x and y
xRange = linspace(-5, 5, 100);
yRange = linspace(-5, 5, 100);

% Create a grid of x and y values
[X, Y] = meshgrid(xRange, yRange);

% Calculate r and theta for each grid point
r = sqrt(X.^2 + Y.^2);
theta = atan2(Y, X);

% Perfectly conducting cylinder
ac(1:51) = besselj((0:50), k0*R)./besselh((0:50), 2, k0*R);
e(1) = 1;
e(2:51) = 2;

Ec = zeros(size(r));
Ec_sc = zeros(size(r));

ll = 1;
s_c = zeros();
R_c = zeros();
Theta_c = zeros();
for ri = 1:numel(r)
    % Calculate field outside of the cylinder
    if r(ri) > R
        Ec_sc(ri) = -sum(e(1:51) .* (-1i).^(0:50) .* ac(1:51) .* besselh(0:50, 2, k0*r(ri)) ...
            .* cos((0:50).*theta(ri)), 2);

        Ec(ri) = Ec_sc(ri) + exp(-1i*X(ri)*k0);

        % Calculate scattering cross section
        if r(ri) > (2*R)^2 
            s_c(ll) = 2*pi*(r(ri)*abs(Ec_sc(ri))^2);
            Theta_c(ll) = theta(ri);
            ll = ll + 1;
        end
    end
end
s_c = s_c .* 4/k0;


% Dielectric cylinder
er = 1.6;
nc = sqrt(er);

for n = 0:50
    if n == 0
        e(n+1) = 1;
    else
        e(n+1) = 2;
    end

    % Calculate parameters field parameters 
    ad(n+1) = ((besselj(n-1, k0*R) - besselj(n+1, k0*R))*besselj(n, nc*k0*R) - ...
        nc*besselj(n, k0*R)*(besselj(n-1, nc*k0*R) - besselj(n+1, nc*k0*R)))/ ...
        ((besselh(n-1, 2, k0*R) - besselh(n+1, 2, k0*R))*besselj(n, nc*k0*R) - ...
        nc*besselh(n, 2, k0*R)*(besselj(n-1, nc*k0*R) - besselj(n+1, nc*k0*R)));

    bd(n+1) = ((besselj(n-1, k0*R) - besselj(n+1, k0*R))*besselh(n, 2, k0*R) - ...
        besselj(n, k0*R)*(besselh(n-1, 2, k0*R) - besselh(n+1, 2, k0*R)))/ ...
        ((besselh(n-1, 2, k0*R) - besselh(n+1, 2, k0*R))*besselj(n, nc*k0*R) - ...
        nc*besselh(n, 2, k0*R)*(besselj(n-1, nc*k0*R) - besselj(n+1, nc*k0*R)));
end

Ed = zeros(size(r));
Ed_sc = zeros(size(r));

mm = 1;
s_d = zeros();
R_d = zeros();
Theta_d = zeros();
for ri = 1:numel(r)
    % Calculate field outside of the cylinder
    if r(ri) > R
        Ed_sc(ri) = -sum(e(1:51) .* (-1i).^(0:50) .* ad(1:51) .* besselh(0:50, 2, k0*r(ri)) ...
            .* cos((0:50).*theta(ri)), 2);

        Ed(ri) = Ed_sc(ri) + exp(-1i*X(ri)*k0);
    else
        % Calculate field inside the cylinder
        Ed_sc(ri) = sum(e(1:51) .* (-1i).^(0:50) .* bd(1:51) .* besselj(0:50, nc*k0*r(ri)) ...
            .* cos((0:50).*theta(ri)), 2);

        Ed(ri) = Ed_sc(ri) + exp(-1i*X(ri)*k0*nc);
    end

    % Calculate scattering cross section
    s_d(mm) = 2*pi*(r(ri)*abs(Ed_sc(ri))^2);
    Theta_d(mm) = theta(ri);
    mm = mm + 1;
end

% Plot Analytical Solution for Perfect Conductor
figure(1);
subplot(1, 2, 1);
pcolor(X, Y, abs(Ec));
title("Total Field")
shading interp;
subplot(1, 2, 2);
pcolor(X, Y, abs(Ec_sc));
title("Scattered Field");
shading interp;
sgtitle("Scattering on perfectly conducting cylinder (Analytical Solution)");

% Plot Bistable Radar crossection for Perfect Conductor
figure(2);
plot(Theta_c*180/pi, s_c, '.');
title("Bistable Radar crossection");
xlabel("Angle (deg)");

% Plot Analytical Solution for Dielectric cylinder
figure(3);
subplot(1, 2, 1);
pcolor(X, Y, abs(Ed));
title("Total Field");
shading interp;
subplot(1, 2, 2);
pcolor(X, Y, abs(Ed_sc));
title("Scattered Field");
shading interp;
sgtitle("Scattering on dielectric cylinder (Analytical Solution)");

% Bistable Radar crossection for Dielectric cylinder
figure(4);
plot(Theta_d*180/pi, s_d, '.');
title("Bistable Radar crossection");
xlabel("Angle (deg)");
