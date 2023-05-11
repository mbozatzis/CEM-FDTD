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
it = 1;
for ti = 0.2:0.2:1
    dt = ti*dx/c0;
    vp = (2./(k*dt)).*asin(c0.*dt.*sin(k*dx/2)./dx);
    plot(k*dx/(2*pi), vp*10^(-8));
    Legend{it}=strcat('|c0*dt/dx| = ', num2str(ti));
    it = it + 1;
end
ylabel("vp (x10^8 m/s)");
xlabel("kdx (pi)")
legend(Legend);

figure(2);
hold on;
theta = 0:0.001:pi/2;
c0 = 3*10^8;
it = 1;
for dxi = 0.02:0.02:0.1
    dx = dxi*lambda0;
    dt = 0.9*dx/(c0*sqrt(2));
    vp = (2./(k*dt)).*asin(c0.*dt.*sqrt(sin(k.*cos(theta).*dx./2).^2+sin(k.*sin(theta).*dx./2).^2)./dx);
    plot(theta*180/pi, vp/c0);
    Legend{it}=strcat('dx = ', num2str(dxi) ,'*lambda0');
    it = it + 1;
end
ylabel("vp/c");
xlabel("theta (deg)")
legend(Legend);
