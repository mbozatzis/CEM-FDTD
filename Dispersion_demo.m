clc
clear

f0 = 10*10^9;
omega0 = 2*pi*f0;
c0 = 3*10^8;

k = omega0/c0;
lambda0 = c0/f0;

dx = (lambda0/100):(lambda0/100):lambda0/2;

figure(1);
hold on;
for ti = 0.1:0.1:1
    dt = ti*dx/c0;
    vp = (2./(k*dt)).*asin(c0.*dt.*sin(k*dx/2)./dx);
    plot(k*dx/pi, vp*10^(-8));
end

figure(2);
theta = 0:0.001:pi/2;
dx = lambda0/20;
dt = 0.7*dx/c0;
vp = (2./(k*dt)).*asin(c0.*dt.*sqrt(sin(k.*cos(theta).*dx./2).^2+sin(k.*sin(theta).*dx./2).^2)./dx);
plot(theta, vp);

% cylinder_options = [3, 1, 3.4, 1.2, 0]; % [x0, R, e_r, sigma_r, y0]
% simulation_options = [10, 10, 12, 10*10^9]; % [X0, Y0, Tn, f0]
% boundary = "PML"; % Type of boundary condition
% boundary_case = "full"; % boundaries only in the left side
% PML_options = [8, 2, 10^(-6)]; % [Npml, power, R]
% 
% [Ez_p, Hx_p, Hy_p] = CylinderScattering(cylinder_options, simulation_options, ...
%     boundary, boundary_case, PML_options);
% 
% for k = 1:size(Ez_p, 3)
%     pcolor(Ez_p(:, :, k));
%     shading interp;
%     colormap default;
%     drawnow;
% end