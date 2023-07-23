clc
clear

% Extract options
Xm_nl = 10;
Ym_nl = 10;
N_x = 100;
N_y = 100;
Tn = 12;
f0 = 2*10^9;

% Useful constants
e0=8.85418781762e-12;
m0=1.256637061436e-06;
c0 = 3*10^8;
lambda0 = c0/f0;

% Simulation Characteristics
X_m = Xm_nl*lambda0;
Y_m = Ym_nl*lambda0;
dx = X_m/N_x;
dy = Y_m/N_y;

Tmax = Tn/f0;
dt = 0.9*dx/(sqrt(2)*c0);

Npml = 8;
pow = 2;
Rpml = 10^(-6);
N_x_new = N_x + 2*Npml;
N_y_new = N_y + 2*Npml;

Ntfsf = 8;
tf_region_start = [Npml+Ntfsf, 1];   
tf_region_end = [N_x+Npml-Ntfsf, N_x_new]; 
fi = 0;

e = e0*ones(N_x_new+1, N_y_new+1);
sigma = zeros(N_x_new+1, N_y_new+1);

Ca = (2*e-dt*sigma)./(2*e+dt*sigma);
Cb = (2*dt./(2*e + dt * sigma))*(1/dx);
Da = dt/(m0*dy);
Db = dt/(m0*dx);

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

% PML Parameters    
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


Cax_pml = exp(1).^(-sigmaEx.*dt/e0);
Cbx_pml = (1-Cax_pml)./(sigmaEx.*dx);
Cay_pml = exp(1).^(-sigmaEy.*dt/e0);
Cby_pml = (1-Cay_pml)./(sigmaEy.*dx);
Da_pml = exp(1).^( -sigmaHz.*dt./m0);
Db_pml = (1-Da_pml)./(sigmaHz.*dx);



for t = 0:dt:Tmax

    % Update Electric Field
    Ex = updateEx(Ex, Hz, Ca, Cb, N_x, N_y, Npml);
    Ey = updateEy(Ey, Hz, Ca, Cb, N_x, N_y, Npml);

    % TFSF Electric (-x)
    Ey(tf_region_start(1), tf_region_start(2):tf_region_end(2)) = ...
            Ey(tf_region_start(1), tf_region_start(2):tf_region_end(2)) + ...
            dt/(e0*dy)*Hz_inc(tf_region_start(1), tf_region_start(2):tf_region_end(2));


    % TFSF Electric (+x)
    Ey(tf_region_end(1), tf_region_start(2):tf_region_end(2)) = ...
            Ey(tf_region_end(1), tf_region_start(2):tf_region_end(2)) - ...
            dt/(e0*dy)*Hz_inc(tf_region_end(1), tf_region_start(2):tf_region_end(2));

    % Update Auxilary field E
    E_aux_prev_max = E_aux(aux_size);
    for i = 2:aux_size
        E_aux(i) = E_aux(i) + dt/(e0*dx)*(H_aux(i-1) - H_aux(i));
    end

    % Auxilary field boundary (Mur 1st)
    E_aux(aux_size+1) = E_aux_prev_max - (dx-c0*dt)/(dx+c0*dt) ...
        *(E_aux(aux_size) - E_aux(aux_size+1));

    if (fi >= 0) && (fi <= pi/2)
        base_i = tf_region_start(1);
        base_j = tf_region_start(2);
    elseif (fi > pi/2) && (fi <= pi)
        base_i = tf_region_end(1);
        base_j = tf_region_start(2);
    elseif (fi > pi) && (fi <= 3*pi/2)
        base_i = tf_region_end(1);
        base_j = tf_region_end(2);
    elseif (fi > 3*pi/2) && (fi <= 2*pi)
        base_i = tf_region_start(1);
        base_j = tf_region_end(2);
    end

    % Calculate incident field Ex from Auxilary field E - xn, xp
    for j = tf_region_start(2):tf_region_end(2)
        d1 = (tf_region_start(1) - base_i)*cos(fi) + (j - base_j)*sin(fi);
        d2 = (tf_region_end(1) - base_i)*cos(fi) + (j - base_j)*sin(fi);
        d1t = d1 - ceil(d1);
        d2t = d2 - ceil(d2);

        Ex_inc(tf_region_start(1), j) = ((1-d1t)*E_aux(2+ceil(d1)) + ...
            d1t*E_aux(2+ceil(d1)+1))*sin(fi);
        Ex_inc(tf_region_end(1), j) = ((1-d2t)*E_aux(2+ceil(d2)) + ...
            d2t*E_aux(2+ceil(d2)+1))*sin(fi);
    end

    % Calculate incident field Ey from Auxilary field E - xn, xp
    for j = tf_region_start(2):tf_region_end(2)
        d1 = (tf_region_start(1) - base_i)*cos(fi) + (j - base_j)*sin(fi);
        d2 = (tf_region_end(1) - base_i)*cos(fi) + (j - base_j)*sin(fi);
        d1t = d1 - ceil(d1);
        d2t = d2 - ceil(d2);

        Ey_inc(tf_region_start(1), j) = ((1-d1t)*E_aux(2+ceil(d1)) + ...
            d1t*E_aux(2+ceil(d1)+1))*(cos(fi));
        Ey_inc(tf_region_end(1), j) = ((1-d2t)*E_aux(2+ceil(d2)) + ...
            d2t*E_aux(2+ceil(d2)+1))*(cos(fi));
    end

    % Calculate incident field Ex from Auxilary field E - yn, yp
    for i = tf_region_start(1):tf_region_end(1)
        d1 = (i - base_i)*cos(fi) + (tf_region_start(2) - base_j)*sin(fi);
        d2 = (i - base_i)*cos(fi) + (tf_region_end(2) - base_j)*sin(fi);
        d1t = d1 - ceil(d1);
        d2t = d2 - ceil(d2);

        Ex_inc(i, tf_region_start(2)) = ((1-d1t)*E_aux(2+ceil(d1)) + ...
            d1t*E_aux(2+ceil(d1)+1))*sin(fi);
        Ex_inc(i, tf_region_end(2)) = ((1-d2t)*E_aux(2+ceil(d2)) + ...
            d2t*E_aux(2+ceil(d2)+1))*sin(fi);
    end

    % Calculate incident field Ey from Auxilary field E - yn, yp
    for i = tf_region_start(1):tf_region_end(1)
        d1 = (i - base_i)*cos(fi) + (tf_region_start(2) - base_j)*sin(fi);
        d2 = (i - base_i)*cos(fi) + (tf_region_end(2) - base_j)*sin(fi);
        d1t = d1 - ceil(d1);
        d2t = d2 - ceil(d2);

        Ey_inc(i, tf_region_start(2)) = ((1-d1t)*E_aux(2+ceil(d1)) + ...
            d1t*E_aux(2+ceil(d1)+1))*(cos(fi));
        Ey_inc(i, tf_region_end(2)) = ((1-d2t)*E_aux(2+ceil(d2)) + ...
            d2t*E_aux(2+ceil(d2)+1))*(cos(fi));
    end


    % PML Electric (-x)
    for i = 2:Npml
        for j = 2:N_y+2*Npml
            Ex(i, j) = Cax_pml(i, j)*Ex(i, j) + Cbx_pml(i, j)*(Hzx_pml(i, j) ...
                + Hzy_pml(i, j) - Hzx_pml(i, j-1) - Hzy_pml(i, j-1));
            Ey(i, j) = Cay_pml(i, j)*Ey(i, j) + Cby_pml(i, j)*(Hzx_pml(i, j) ...
                + Hzy_pml(i, j) - Hzx_pml(i-1, j) - Hzy_pml(i-1, j));
        end
    end

    %PML Electric (+x)
    for i = Npml+N_x+1:N_x+2*Npml-1
        for j = 2:N_y+2*Npml
            Ex(i, j) = Cax_pml(i, j)*Ex(i, j) + Cbx_pml(i, j)*(Hzx_pml(i, j) ...
                + Hzy_pml(i, j) - Hzx_pml(i, j-1) - Hzy_pml(i, j-1));
            if i ~= Npml+N_x+1
                Ey(i, j) = Cay_pml(i, j)*Ey(i, j) + Cby_pml(i, j)*(Hzx_pml(i, j) ...
                    + Hzy_pml(i, j) - Hzx_pml(i-1, j) - Hzy_pml(i-1, j));
            else
                Ey(i, j) = Cay_pml(i, j)*Ey(i, j) + Cby_pml(i, j)*(Hzx_pml(i, j) ...
                    + Hzy_pml(i, j) - Hz(Npml+N_x, j));
            end
        end
    end

    % PBC Electric
    for i = 1:N_x_new
        Ex(i, N_x_new) = Ex(i, N_x_new) + Cb(i, j)*(Hz(i, 1)-Hz(i, N_x_new));
        % Ex(i, Npml) = Ex(i, Npml) + Cb(i, Npml)*(Hz(i, Npml+1)-Hz(i, Npml));
        % Ex(i, Npml+N_x+1) = Ex(i, Npml);
        % 
        % Ey(i, Npml) = Ey(i, Npml) + Cb(i, Npml)*(Hz(i, Npml)-Hz(i+1, Npml));
        % Ey(i, Npml+N_x+1) = Ey(i, Npml);

        Ex(i, N_y_new) = Ex(i, 1);
        Ey(i, N_y_new) = Ey(i, 1);
    end


%     %PML Electric (-y)
%     for i = Npml+1:N_x+Npml
%         for j = 2:Npml
%             Ex(i, j) = Cax_pml(i, j)*Ex(i, j) + Cbx_pml(i, j)*(Hzx_pml(i, j) ...
%                 + Hzy_pml(i, j) - Hzx_pml(i, j-1) - Hzy_pml(i, j-1));
%             Ey(i, j) = Cay_pml(i, j)*Ey(i, j) + Cby_pml(i, j)*(Hzx_pml(i, j) ...
%                 + Hzy_pml(i, j) - Hzx_pml(i-1, j) - Hzy_pml(i-1, j));
%         end
%     end
% 
%     %PML Electric (+y)
%     for i = Npml+1:N_y+Npml
%         for j = Npml+N_x+1:N_x+2*Npml-1
%             if j ~= Npml+N_x+1
%                 Ex(i, j) = Cax_pml(i, j)*Ex(i, j) + Cbx_pml(i, j)*(Hzx_pml(i, j) ...
%                     + Hzy_pml(i, j) - Hzx_pml(i, j-1) - Hzy_pml(i, j-1));
%             else
%                 Ex(i, j) = Cay_pml(i, j)*Ex(i, j) + Cby_pml(i, j)*(Hzx_pml(i, j) ...
%                 + Hzy_pml(i, j) - Hz(i, Npml+N_y));
%             end
%             Ey(i, j) = Cay_pml(i, j)*Ey(i, j) + Cby_pml(i, j)*(Hzx_pml(i, j) ...
%                 + Hzy_pml(i, j) - Hzx_pml(i-1, j) - Hzy_pml(i-1, j));
%         end
%     end

    % Update Magnetic Field
    Hz = updateHz(Hz, Ex, Ey, Da, Db, N_x, N_y, Npml);

    % TFSF Magnetic (-x)
    Hz(tf_region_start(1), tf_region_start(2):tf_region_end(2)) = ...
        Hz(tf_region_start(1), tf_region_start(2):tf_region_end(2)) + ...
        dt/(m0*dy)*Ey_inc(tf_region_start(1), tf_region_start(2):tf_region_end(2));


    % TFSF Magnetic (+x)
    Hz(tf_region_end(1), tf_region_start(2):tf_region_end(2)) = ...
        Hz(tf_region_end(1), tf_region_start(2):tf_region_end(2)) - ...
        dt/(m0*dy)*Ey_inc(tf_region_end(1), tf_region_start(2):tf_region_end(2));

    % Update Auxilary field H
    for i = 1:aux_size
        H_aux(i) = H_aux(i) + dt/(m0*dx)*(E_aux(i) - E_aux(i+1));
    end

    H_aux(1) = sin(2*pi*f0*t);

    % Calculate r components for each angle
    if (fi >= 0) && (fi <= pi/2)
        base_i = tf_region_start(1);
        base_j = tf_region_start(2);
    elseif (fi > pi/2) && (fi <= pi)
        base_i = tf_region_end(1)-1/2;
        base_j = tf_region_start(2);
    elseif (fi > pi) && (fi <= 3*pi/2)
        base_i = tf_region_end(1)-1/2;
        base_j = tf_region_end(2)-1/2;
    elseif (fi > 3*pi/2) && (fi <= 2*pi)
        base_i = tf_region_start(1);
        base_j = tf_region_end(2)-1/2;
    end

    % Calculate incident field Hz from Auxilary field H - xn, xp
    for j = tf_region_start(2):tf_region_end(2)
        d1 = (tf_region_start(1) - base_i)*cos(fi) + (j - base_j)*sin(fi);
        d2 = (tf_region_end(1) - base_i)*cos(fi) + (j - base_j)*sin(fi);
        d1t = d1 - ceil(d1);
        d2t = d2 - ceil(d2);

        Hz_inc(tf_region_start(1), j) = (1-d1t)*H_aux(2+ceil(d1)) + d1t*H_aux(2+ceil(d1)+1);
        Hz_inc(tf_region_end(1), j) = (1-d2t)*H_aux(2+ceil(d2)) + d2t*H_aux(2+ceil(d2)+1);
    end

    % Calculate incident field Hz from Auxilary field H - yn, yp
    for i = tf_region_start(1):tf_region_end(1)
        d1 = (i - base_i)*cos(fi) + (tf_region_start(2) - base_j)*sin(fi);
        d2 = (i - base_i)*cos(fi) + (tf_region_end(2) - base_j)*sin(fi);
        d1t = d1 - ceil(d1);
        d2t = d2 - ceil(d2);

        Hz_inc(i, tf_region_start(2)) = (1-d1t)*H_aux(2+ceil(d1)) + d1t*H_aux(2+ceil(d1)+1);
        Hz_inc(i, tf_region_end(2)) = (1-d2t)*H_aux(2+ceil(d2)) + d2t*H_aux(2+ceil(d2)+1);
    end

    % Update source
    % Hz(floor(N_x_new/2), floor(N_y_new/2)) = sin(2*pi*f0*t);

    % PML Magnetic (-x)
    for i = 2:Npml
        for j = 1:N_y+2*Npml-1
            Hzx_pml(i, j) = Da_pml(i, j)*Hzx_pml(i, j) + Db_pml(i, j)*(Ey(i+1, j) ...
                - Ey(i, j));
            Hzy_pml(i, j) = Da_pml(i, j)*Hzy_pml(i, j) + Db_pml(i, j)*(Ex(i, j+1) ...
                - Ex(i, j));
            Hz(i, j) = Hzx_pml(i, j) + Hzy_pml(i, j);
        end
    end

    %PML Magnetic (+x)
    for i = Npml+N_x+1:N_x+2*Npml-1
        for j = 1:N_y+2*Npml-1
            Hzx_pml(i, j) = Da_pml(i, j)*Hzx_pml(i, j) + Db_pml(i, j)*(Ey(i+1, j) ...
                - Ey(i, j));
            Hzy_pml(i, j) = Da_pml(i, j)*Hzy_pml(i, j) + Db_pml(i, j)*(Ex(i, j+1) ...
                - Ex(i, j));
            Hz(i, j) = Hzx_pml(i, j) + Hzy_pml(i, j);
        end
    end

    % PBC Magnetic
    for i = 2:N_x_new
        Hz(i, 1) = Hz(i, 1) + Db*(Ex(i, 1)-Ex(i, N_y_new) ...
            + Ey(i-1,1) - Ey(i, 1));

        Hz(i, N_y_new) = Hz(i, 1);
    end

%     %PML Magnetic (-y)
%     for i = Npml+1:N_x+Npml
%         for j = 2:Npml
%             Hzx_pml(i, j) = Da_pml(i, j)*Hzx_pml(i, j) + Db_pml(i, j)*(Ey(i+1, j) ...
%                 - Ey(i, j));
%             Hzy_pml(i, j) = Da_pml(i, j)*Hzy_pml(i, j) + Db_pml(i, j)*(Ex(i, j+1) ...
%                 - Ex(i, j));
%             Hz(i, j) = Hzx_pml(i, j) + Hzy_pml(i, j);
%         end
%     end 

%     %PML Magnetic (+y)
%     for i = Npml+1:N_y+Npml
%         for j = Npml+N_x+1:N_x+2*Npml-1
%             Hzx_pml(i, j) = Da_pml(i, j)*Hzx_pml(i, j) + Db_pml(i, j)*(Ey(i+1, j) ...
%                 - Ey(i, j));
%             Hzy_pml(i, j) = Da_pml(i, j)*Hzy_pml(i, j) + Db_pml(i, j)*(Ex(i, j+1) ...
%                 - Ex(i, j));
%             Hz(i, j) = Hzx_pml(i, j) + Hzy_pml(i, j);
%         end
%     end


    pcolor(Hz(Npml+1:Npml+N_y+1,1:2*Npml+N_y));
    shading interp;
    drawnow;
end

