clc
clear

addpath("TE_wave/PML/");
addpath("TE_wave/TFSF/");

% Define options
Xm_nl = 10;
Ym_nl = 10;
N_x = 100;
N_y = 100;
Tn = 15;
f0 = 5*10^9;

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
tf_region_start = [Npml+Ntfsf, Npml+Ntfsf];   
tf_region_end = [N_x+Npml-Ntfsf, N_y+Npml-Ntfsf]; 
fi = pi/4;

% Medium Characteristics
e = e0*ones(N_x_new+1, N_y_new+1);
sigma = zeros(N_x_new+1, N_y_new+1);



% Field updating parameters
Ca = (2*e-dt*sigma)./(2*e+dt*sigma);
Cb = (2*dt./(2*e + dt * sigma))*(1/dx);
Da = dt/(m0*dy);
Db = dt/(m0*dx);

% Initialize fields
Hz = zeros(N_x_new+1, N_y_new+1);
Ex = zeros(N_x_new, N_y_new+1);
Ey = zeros(N_x_new+1, N_y_new);

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
    Ex = updateEx(Ex, Hz, Ca, Cb, N_x, N_y, Npml, 1);
    Ey = updateEy(Ey, Hz, Ca, Cb, N_x, N_y, Npml, 1);

    % Ez Correction on the TFSF boundaries
    [Ex, Ey] = updateTFSFboundE(Ex, Ey, Hz_inc, tf_region_start, tf_region_end, dx, dy, dt);

    % Update Auxilary field E
    E_aux_prev_max = E_aux(aux_size);
    for i = 2:aux_size
        E_aux(i) = E_aux(i) + dt/(e0*dx)*(H_aux(i-1) - H_aux(i));
    end

    % Auxilary field boundary (Mur 1st)
    E_aux(aux_size+1) = E_aux_prev_max - (dx-c0*dt)/(dx+c0*dt) ...
        *(E_aux(aux_size) - E_aux(aux_size+1));

    % Calculate incident field Ez from Auxilary field E
    [Ex_inc, Ey_inc] = IncFromAuxE(E_aux, Ex_inc, Ey_inc, fi, tf_region_start, ...
        tf_region_end);
    

    % PML Electric (-x)
    [Ex, Ey] = updatePMLxnE(Ex, Ey, Hzx_pml, Hzy_pml, Cax_pml, ...
    Cay_pml, Cbx_pml, Cby_pml, Npml, N_x, N_y);

    % PML Electric (+x)
    [Ex, Ey] = updatePMLxpE(Ex, Ey, Hzx_pml, Hzy_pml, Hz, Cax_pml, ...
    Cay_pml, Cbx_pml, Cby_pml, Npml, N_x, N_y);

    % PML Electric (-y)
    [Ex, Ey] = updatePMLynE(Ex, Ey, Hzx_pml, Hzy_pml, Cax_pml, ...
    Cay_pml, Cbx_pml, Cby_pml, Npml, N_x, N_y);

    % PML Electric (+y)
    [Ex, Ey] = updatePMLypE(Ex, Ey, Hzx_pml, Hzy_pml, Hz, Cax_pml, ...
    Cay_pml, Cbx_pml, Cby_pml, Npml, N_x, N_y);


    % Update Magnetic Field
    Hz = updateHz(Hz, Ex, Ey, Da, Db, N_x, N_y, Npml, 1);

    % Update Auxilary field H
    for i = 1:aux_size
        H_aux(i) = H_aux(i) + dt/(m0*dx)*(E_aux(i) - E_aux(i+1));
    end

    % Update source on the Auxilary Field
    H_aux(1) = sin(2*pi*f0*t);

    % H Correction on the TFSF boundaries
    Hz = updateTFSFboundH(Hz, Ex_inc, Ey_inc, tf_region_start, ...
        tf_region_end, dx, dy, dt);

    % Calculate incident field Hx and Hy from Auxilary field H
    Hz_inc = IncFromAuxH(H_aux, Hz_inc, fi, tf_region_start, tf_region_end);

    % PML Magnetic (-x)
    [Hzx_pml, Hzy_pml, Hz] = updatePMLxnH(Hzx_pml, Hzy_pml, ...
    Hz, Ex, Ey, Da_pml, Db_pml, Npml, N_x, N_y);

    % PML Magnetic (+x)
    [Hzx_pml, Hzy_pml, Hz] = updatePMLxpH(Hzx_pml, Hzy_pml, ...
    Hz, Ex, Ey, Da_pml, Db_pml, Npml, N_x, N_y);

    % PML Magnetic (-y)
    [Hzx_pml, Hzy_pml, Hz] = updatePMLynH(Hzx_pml, Hzy_pml, ...
    Hz, Ex, Ey, Da_pml, Db_pml, Npml, N_x, N_y);

    % PML Magnetic (+y)
    [Hzx_pml, Hzy_pml, Hz] = updatePMLypH(Hzx_pml, Hzy_pml, ...
    Hz, Ex, Ey, Da_pml, Db_pml, Npml, N_x, N_y);


    pcolor(real(Hz(Npml+1:Npml+N_y+1, Npml+1:Npml+N_y+1)));
    shading interp;
    drawnow;

end