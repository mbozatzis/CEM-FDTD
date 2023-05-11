clc
clear

Npml = 16;
Rpml = 10^(-6);
pow = 2;
boundary_case = "full";

f0 = 10*10^9;
Xm_nl = 10;
Ym_nl = 10;
N_x = 100;
N_y = 100;
Tn = 12;

e0=8.85418781762e-12;
m0=1.256637061436e-06;
c0 = 3*10^8;
lambda0 = c0/f0;

X_m = Xm_nl*lambda0;
Y_m = Ym_nl*lambda0;

dx = X_m/N_x;
dy = Y_m/N_y;

Tmax = Tn/f0;
dt = 0.9*dx/(sqrt(2)*c0);

N_x_new = N_x + 2*Npml;
N_y_new = N_y + 2*Npml;


e = e0*ones(N_x_new+1, N_y_new+1);
sigma = zeros(N_x_new+1, N_y_new+1);



% TM: Ez, Hx, Hy
Ez = zeros(N_x_new+1, N_y_new+1);
Hx = zeros(N_x_new+1, N_y_new);
Hy = zeros(N_x_new, N_y_new+1);

Ezx_pml = zeros(N_x_new+1, N_y_new+1);
Ezy_pml = zeros(N_x_new+1, N_y_new+1);

% Constants for field updating
Ca = (2*e-dt*sigma)./(2*e+dt*sigma);
Cb = (2*dt./(2*e + dt * sigma))*(1/dx);
Da = dt/(m0*dy);
Db = dt/(m0*dx);


% se= -e0*c0*log(Rpml)/(2^(pow+2)*dx*Npml^(pow+1));
se = - (pow+1)*e0*c0*log(Rpml)/(2*dx*Npml);
sh=(m0/e0)*se;
sigmaE = zeros(N_x_new, N_y_new);
sigmaHx = zeros(N_x_new, N_y_new);
sigmaHy = zeros(N_x_new, N_y_new);

% X-negative (xn)
for i = 1:Npml
%     sigmaE(Npml-i+1, :)= sigmaE(Npml-i+1, :) + se*((2*i+1)^(pow+1)-(2*i-1)^(pow+1));
%     sigmaHx(Npml-i+1, :)= sigmaHx(Npml-i+1,:) + sh*((2*i+1)^(pow+1)-(2*i-1)^(pow+1));
%     sigmaHy(Npml-i+1, :)= sigmaHy(Npml-i+1,:) + sh*((2*(i+0.5)+1)^(pow+1)-(2*(i+0.5)-1)^(pow+1));
    sigmaE(Npml-i+1, :)= sigmaE(Npml-i+1, :) + se*(i/Npml)^(pow);
    sigmaHx(Npml-i+1, :)= sigmaHx(Npml-i+1,:) + sh*(i/Npml)^(pow);
    sigmaHy(Npml-i+1, :)= sigmaHy(Npml-i+1,:) + sh*((i+0.5)/Npml)^(pow);
end

% X-positive (xp)
for i = 1:Npml
%     sigmaE(Npml+N_x+i, :)= sigmaE(Npml+N_x+i, :) + se*((2*i+1)^(pow+1)-(2*i-1)^(pow+1));
%     sigmaHx(Npml+N_x+i, :)= sigmaHx(Npml+N_x+i,:) + sh*((2*i+1)^(pow+1)-(2*i-1)^(pow+1));
%     sigmaHy(Npml+N_x+i, :)= sigmaHy(Npml+N_x+i,:) + sh*((2*(i+0.5)+1)^(pow+1)-(2*(i+0.5)-1)^(pow+1));
    sigmaE(Npml+N_x+i, :)= sigmaE(Npml+N_x+i, :) + se*(i/Npml)^(pow);
    sigmaHx(Npml+N_x+i, :)= sigmaHx(Npml+N_x+i,:) + sh*(i/Npml)^(pow);
    sigmaHy(Npml+N_x+i, :)= sigmaHy(Npml+N_x+i,:) + sh*((i+0.5)/Npml)^(pow);
end 

% Y-negative (yn)
for j = 1:Npml
%     sigmaE(:, Npml-j+1)= sigmaE(:, Npml-j+1) + se*((2*j+1)^(pow+1)-(2*j-1)^(pow+1));
%     sigmaHx(:, Npml-j+1)= sigmaHx(:, Npml-j+1) + sh*((2*j+1)^(pow+1)-(2*j-1)^(pow+1));
%     sigmaHy(:, Npml-j+1)= sigmaHy(:, Npml-j+1) + sh*((2*(j+0.5)+1)^(pow+j)-(2*(j+0.5)-1)^(pow+1));
    sigmaE(:, Npml-j+1)= sigmaE(:, Npml-j+1) + se*(j/Npml)^(pow);
    sigmaHx(:, Npml-j+1)= sigmaHx(:, Npml-j+1) + sh*(j/Npml)^(pow);
    sigmaHy(:, Npml-j+1)= sigmaHy(:, Npml-j+1) + sh*((j+0.5)/Npml)^(pow);
end

% Y-positive (yp)
for j = 1:Npml
%     sigmaE(:, Npml+N_y+j)= sigmaE(:, Npml+N_y+j) + se*((2*j+1)^(pow+1)-(2*j-1)^(pow+1));
%     sigmaHx(:, Npml+N_y+j)= sigmaHx(:, Npml+N_y+j) + sh*((2*j+1)^(pow+1)-(2*j-1)^(pow+1));
%     sigmaHy(:, Npml+N_y+j)= sigmaHy(:, Npml+N_y+j) + sh*((2*(j+0.5)+1)^(pow+j)-(2*(j+0.5)-1)^(pow+1));
    sigmaE(:, Npml+N_y+j)= sigmaE(:, Npml+N_y+j) + se*(j/Npml)^(pow);
    sigmaHx(:, Npml+N_y+j)= sigmaHx(:, Npml+N_y+j) + sh*(j/Npml)^(pow);
    sigmaHy(:, Npml+N_y+j)= sigmaHy(:, Npml+N_y+j) + sh*((j+0.5)/Npml)^(pow);
end 

Ca_pml = exp(1).^(-sigmaE.*dt/e0);
Cb_pml = (1-Ca_pml)./(sigmaE.*dx);
Dax_pml = exp(1).^( -sigmaHx.*dt./m0);
Day_pml = exp(1).^( -sigmaHy.*dt./m0);
Dbx_pml = (1-Dax_pml)./(sigmaHx.*dx);
Dby_pml = (1-Day_pml)./(sigmaHy.*dx);

l = 1;
for t = 0:dt:Tmax            
    % Update Electric field
    Ez = updateEz(Ez, Hx, Hy, Ca, Cb, N_x, N_y, 1, Npml);

    % Update source
    Ez(floor(N_x_new/2), floor(N_y_new/2)) = sin(2*pi*f0*t);       

    % Update PML (-x) Electric
    [Ezx_pml, Ezy_pml] = updatePMLxnE(Ezx_pml, Ezy_pml, Hx, Hy, Ca_pml, Cb_pml, Npml, N_x, N_y);   

    % Update PML (+x) Electric
    [Ezx_pml, Ezy_pml] = updatePMLxpE(Ezx_pml, Ezy_pml, Hx, Hy, Ca_pml, Cb_pml, Npml, N_x, N_y);

    % Update PML (-y) Electric 
    [Ezx_pml, Ezy_pml] = updatePMLynE(Ezx_pml, Ezy_pml, Hx, Hy, Ca_pml, Cb_pml, Npml, N_x, N_y);

    % Update PML (+y) Electric
    [Ezx_pml, Ezy_pml] = updatePMLypE(Ezx_pml, Ezy_pml, Hx, Hy, Ca_pml, Cb_pml, Npml, N_x, N_y);

    % Update Magnetic Field
    Hx = updateHx(Hx, Ez, Da, N_x, N_y, 1, Npml);
    Hy = updateHy(Hy, Ez, Db, N_x, N_y, 1, Npml);
    
    % Update PML (-x) Magnetic
    [Hx, Hy] = updatePMLxnH(Hx, Hy, Ez, Ezx_pml, Ezy_pml, Dax_pml, Dbx_pml, Day_pml, ...
        Dby_pml, Npml, N_x, N_y);

    % Update PML (+x) Magnetic
    [Hx, Hy] = updatePMLxpH(Hx, Hy, Ez, Ezx_pml, Ezy_pml, Dax_pml, Dbx_pml, Day_pml, ...
        Dby_pml, Npml, N_x, N_y);

    % Update PML (-y) Magnetic
    [Hx, Hy] = updatePMLynH(Hx, Hy, Ez, Ezx_pml, Ezy_pml, Dax_pml, Dbx_pml, Day_pml, ...
        Dby_pml, Npml, N_x, N_y);

    % Update PML (+y) Magnetic
    [Hx, Hy] = updatePMLypH(Hx, Hy, Ez, Ezx_pml, Ezy_pml, Dax_pml, Dbx_pml, Day_pml, ...
        Dby_pml, Npml, N_x, N_y);

    eez(:, :, l) = Ez(Npml+1:Npml+N_y+1,Npml+1:Npml+N_y+1);
    l = l + 1;
end




cylinder_options = [0, 0, 1, 0, 0];
PML_options = [8, 2, 10^(-6)]; 
% Reference Field simulation
boundary = "No-boundary";
simulation_options = [20, 20, 200, 200, 12, 10*10^9];
[Ez_ref, Hx_ref, Hy_ref] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

% Add testpoint at the middle of the bottom boundary
test_point = [1, floor(size(eez, 1)/2)]; 

% Add the tespoint at the relevant position in the reference grid
ref_test_point = [floor(size(Ez_ref, 1)/2), floor(size(Ez_ref, 2)/2) - floor(size(eez, 2)/2) + 1];


Et_p = reshape(eez(test_point(1), test_point(2), :), [1, size(eez, 3)]);
Et_ref = reshape(Ez_ref(ref_test_point(2), ref_test_point(1), :), [1, size(Ez_ref, 3)]);

figure(1);
pcolor(eez(:, :, end));
shading interp;
hold on;
xline(test_point(1));
yline(test_point(2));

figure(2);
pcolor(Ez_ref(:, :, end));
shading interp;
hold on;
xline(ref_test_point(2));
yline(ref_test_point(1));

ref = max(Et_ref);
figure(3);
sgtitle("Relative error");
plot(abs(Et_ref-Et_p)/ref);


