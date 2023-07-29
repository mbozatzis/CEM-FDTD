clc
clear

addpath("TE_wave/PML/");

% Define options
Xm_nl = 10;
Ym_nl = 10;
N_x = 100;
N_y = 100;
Tn = 7;
f0 = 0.5*10^12;

% Useful constants
e0=8.85418781762e-12;
m0=1.256637061436e-06;
c0 = 3*10^8;
lambda0 = c0/f0;
k0 = 2*pi*f0/c0;

% Simulation Characteristics
X_m = Xm_nl*lambda0;
Y_m = Ym_nl*lambda0;
dx = X_m/N_x;
dy = Y_m/N_y;

Tmax = Tn/f0;
dt = 0.9*dx/(sqrt(2)*c0);

% PML Characteristics
Npml = 8;
pow = 2;
Rpml = 10^(-6);
N_x_new = N_x + 2*Npml;
N_y_new = N_y + 2*Npml;

% TFSF Characteristics
Ntfsf = 8;
tf_region_start = [Npml+Ntfsf, 1];   
tf_region_end = [N_x+Npml-Ntfsf, N_x_new]; 
fi = pi;

% Medium Characteristics
e = e0*ones(N_x_new+1, N_y_new+1);
sigma = zeros(N_x_new+1, N_y_new+1);
e1 = e0;
e2 = e0;

% e1 = e0;
% e2 = e0;
% for i = 1:N_x_new+1
%     for j = 1:N_y_new+1
%         if i>= 58
%             e(i, j) = e(i, j);
%         end
%     end
% end

% Graphene Characteristics
A_gr = 0.1166; %0.02344;
tau_gr = 3.3151e-12; %6.648e-13;
pos_gr = Npml + N_x/2; 


% Field updating parameters
Ca = (2*e-dt*sigma)./(2*e+dt*sigma);
Cb = (2*dt./(2*e + dt * sigma))*(1/dx);
Da = dt/(m0*dy);
Db = dt/(m0*dx);
Cb_gr = 2*dt/((e1 + e2)*dx);

% Initialize fields
Hz = zeros(N_x_new+1, N_y_new+1);
Ex = zeros(N_x_new, N_y_new+1);
Ey = zeros(N_x_new+1, N_y_new);
Jy = zeros(N_x_new+1, N_y_new);

Hzx_pml = zeros(N_x_new+1, N_y_new+1);
Hzy_pml = zeros(N_x_new+1, N_y_new+1);

Hz_inc = zeros(N_x_new+1, N_y_new+1);
Ex_inc = zeros(N_x_new+1, N_y_new+1);
Ey_inc = zeros(N_x_new+1, N_y_new+1);

% Define Auxilary Field for TFSF
aux_size = ceil(N_x*sqrt(2)+4);
E_aux = zeros(aux_size + 1, 1);
H_aux = zeros(aux_size, 1);

% PML Conductivities    
se = - (pow+1)*e0*c0*log(Rpml)/(2*dx*Npml);
sh = (m0/e0)*se;
sigmaEx = zeros(N_x_new, N_y_new);
sigmaEy = zeros(N_x_new, N_y_new);
sigmaHz = zeros(N_x_new, N_y_new);

% X-negative (-x)
for i = 1:Npml
    sigmaEx(Npml-i+1, :)= sigmaEx(Npml-i+1, :) + se*((i+0.5)/Npml)^(pow);
    sigmaEy(Npml-i+1, :)= sigmaEy(Npml-i+1, :) + se*(i/Npml)^(pow);
    sigmaHz(Npml-i+1, :)= sigmaHz(Npml-i+1,:) + sh*((i+0.5)/Npml)^(pow);
end

% X-positive (+x)
for i = 1:Npml
    sigmaEx(Npml+N_x+i, :)= sigmaEx(Npml+N_x+i, :) + se*((i+0.5)/Npml)^(pow);
    sigmaEy(Npml+N_x+i, :)= sigmaEy(Npml+N_x+i, :) + se*(i/Npml)^(pow);
    sigmaHz(Npml+N_x+i, :)= sigmaHz(Npml+N_x+i,:) + sh*((i+0.5)/Npml)^(pow);
end 

% Y-negative (-y)
for j = 1:Npml
    sigmaEx(:, Npml-j+1)= sigmaEx(:, Npml-j+1) + se*(j/Npml)^(pow);
    sigmaEy(:, Npml-j+1)= sigmaEy(:, Npml-j+1) + se*((j+0.5)/Npml)^(pow);
    sigmaHz(:, Npml-j+1)= sigmaHz(:, Npml-j+1) + sh*((j+0.5)/Npml)^(pow);
end

% Y-positive (+y)
for j = 1:Npml
    sigmaEx(:, Npml+N_y+j)= sigmaEx(:, Npml+N_y+j) + se*(j/Npml)^(pow);
    sigmaEy(:, Npml+N_y+j)= sigmaEy(:, Npml+N_y+j) + se*((j+0.5)/Npml)^(pow);
    sigmaHz(:, Npml+N_y+j)= sigmaHz(:, Npml+N_y+j) + sh*((j+0.5)/Npml)^(pow);
end 


% PML updating parameters
Cax_pml = exp(1).^(-sigmaEx.*dt/e0);
Cbx_pml = (1-Cax_pml)./(sigmaEx.*dx);
Cay_pml = exp(1).^(-sigmaEy.*dt/e0);
Cby_pml = (1-Cay_pml)./(sigmaEy.*dx);
Da_pml = exp(1).^( -sigmaHz.*dt./m0);
Db_pml = (1-Da_pml)./(sigmaHz.*dx);


for t = 0:dt:Tmax

    % Update Electric Field
    Ex = updateEx(Ex, Hz, Ca, Cb, N_x, N_y, Npml);
    Ey = updateGrapheneEy(Ey, Hz, Ca, Cb, N_x, N_y, Npml, pos_gr, Cb_gr, Jy);

    % Update Graphene ADE
    Jy = updateGrapheneADE(Jy, Ey, tau_gr, A_gr, dt);
    
    % Update source
    Ey(floor(N_x_new/2+1), floor(N_y_new/2)) = sin(2*pi*f0*t);

    % PML Electric (-x)
    [Ex, Ey] = updatePMLxnE(Ex, Ey, Hzx_pml, Hzy_pml, Cax_pml, ...
    Cay_pml, Cbx_pml, Cby_pml, Npml, N_x, N_y);

    %PML Electric (+x)
    [Ex, Ey] = updatePMLxpE(Ex, Ey, Hzx_pml, Hzy_pml, Hz, Cax_pml, ...
    Cay_pml, Cbx_pml, Cby_pml, Npml, N_x, N_y);

    % PBC Electric
    for i = 1:N_x_new
        Ex(i, N_x_new) = Ex(i, N_x_new) + Cb(i, j)*(Hz(i, 1)-Hz(i, N_x_new));
        % Ex(i, Npml) = Ex(i, Npml) + Cb(i, Npml)*(Hz(i, Npml+1)-Hz(i, Npml));
        % Ex(i, Npml+N_x+1) = Ex(i, Npml);
        % 
        % Ey(i, Npml) = Ey(i, Npml) + Cb(i, Npml)*(Hz(i, Npml)-Hz(i+1, Npml));
        % Ey(i, Npml+N_x+1) = Ey(i, Npml);

        Ex(i, 1) = Ex(i, N_x_new)*exp(1i*k0*N_y_new);
        Ey(i, 1) = Ey(i, N_x_new)*exp(1i*k0*N_y_new);
    end


    % Update Magnetic Field
    Hz = updateHz(Hz, Ex, Ey, Da, Db, N_x, N_y, Npml);

    % PML Magnetic (-x)
    [Hzx_pml, Hzy_pml, Hz] = updatePMLxnH(Hzx_pml, Hzy_pml, ...
    Hz, Ex, Ey, Da_pml, Db_pml, Npml, N_x, N_y);

    %PML Magnetic (+x)
    [Hzx_pml, Hzy_pml, Hz] = updatePMLxpH(Hzx_pml, Hzy_pml, ...
    Hz, Ex, Ey, Da_pml, Db_pml, Npml, N_x, N_y);

    % PBC Magnetic
    for i = 2:N_x_new
        Hz(i, 1) = Hz(i, 1) + Db*(Ex(i, 1)-Ex(i, N_y_new) ...
            + Ey(i-1,1) - Ey(i, 1));

        Hz(i, N_y_new) = Hz(i, 1)*exp(-1i*k0*N_y_new);
    end


    pcolor(real(Hz(Npml+1:Npml+N_y+1,1:2*Npml+N_y)));
    shading interp;
    drawnow;

end

