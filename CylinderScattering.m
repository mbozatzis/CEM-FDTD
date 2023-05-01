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
e_r = 1;%3.4;
sigma_r = 0;%1.2;

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
Cb = (2*dt./(2*e + dt * sigma))*(1/dx);
Da = dt/(m0*dy);
Db = dt/(m0*dx);

% % Without boundary conditions (zero)
% for t = 0:dt:Tmax
%     for i = 2:N_x
%         for j = 2:N_y
%             Ez(i, j) = Ca(i, j)*Ez(i, j) + Cb(i, j)*(Hy(i, j)-Hy(i-1, j) ...
%                 + Hx(i,j-1) - Hx(i, j));
%         end
%     end
% 
% 
%     Ez(N_x/2, N_y/2) = sin(2*pi*f0*t);
% 
%     for i = 2:N_x
%         for j = 1:N_y
%             Hx(i, j) = Hx(i, j) - Da*(Ez(i, j+1) - Ez(i, j));
%         end
%     end
% 
%     for i = 1:N_x
%         for j = 2:N_y
%             Hy(i, j) = Hy(i, j) + Db*(Ez(i+1, j) - Ez(i, j));
%         end
%     end
%     
%     surf(Ez);
%     colormap default;
%     drawnow;
% 
% end


% % With first order Mur boundary conditions
% for t = 0:dt:Tmax
%     for i = 1:N_x+1
%         E2yprev = Ez(i, 2);
%         ENyprev = Ez(i, N_y);
%         for j = 1:N_y+1
%             E2xprev = Ez(2, j);
%             ENxprev = Ez(N_x, j);
%             if (i ~= 1) && (j ~= 1) && (i ~= N_x+1) && (j ~= N_y+1)
%                 Ez(i, j) = Ca(i, j)*Ez(i, j) + Cb(i, j)*(Hy(i, j)-Hy(i-1, j) ...
%                     + Hx(i,j-1) - Hx(i, j));
%             end
%             Ez(1, j) = E2xprev + ((c0*dt-dx)/(c0*dt+dx))*(Ez(2, j) - Ez(1, j));
%             Ez(N_x+1, j) = ENxprev + ((c0*dt-dx)/(c0*dt+dx))*(Ez(N_x, j) - Ez(N_x+1, j));
%         end
%         Ez(i, 1) = E2yprev + ((c0*dt-dx)/(c0*dt+dx))*(Ez(i, 2) - Ez(i, 1));
%         Ez(i, N_y+1) = ENyprev + ((c0*dt-dx)/(c0*dt+dx))*(Ez(i, N_y) - Ez(i, N_y+1));
%     end
% 
% 
%     Ez(N_x/2, N_y/2) = sin(2*pi*f0*t);
% 
%     for i = 1:N_x+1
%         for j = 1:N_y
%             Hx(i, j) = Hx(i, j) - Da*(Ez(i, j+1) - Ez(i, j));
%         end
%     end
% 
%     for i = 1:N_x
%         for j = 1:N_y+1
%             Hy(i, j) = Hy(i, j) + Db*(Ez(i+1, j) - Ez(i, j));
%         end
%     end
%     
%     surf(Ez);
%     colormap default;
%     drawnow;
% 
% end


% % With second order Mur boundary conditions
% 
% Ez1_t1_x = Ez(1,:);
% Ez2_t1_x = Ez(2, :);
% Ezmax_t1_x = Ez(N_x+1,:);
% Ezmax2_t1_x = Ez(N_x, :);
% Ez1_t1_y = Ez(:,1);
% Ez2_t1_y = Ez(:, 2);
% Ezmax_t1_y = Ez(:, N_y+1);
% Ezmax2_t1_y = Ez(:, N_y);
% 
% for t = 0:dt:Tmax
%     Ez1_t2_x = Ez1_t1_x;
%     Ez1_t1_x = Ez(1,:);
%     Ez2_t2_x = Ez2_t1_x;
%     Ez2_t1_x = Ez(2, :);
% 
%     Ezmax_t2_x = Ezmax_t1_x;
%     Ezmax_t1_x = Ez(N_x+1,:);
%     Ezmax2_t2_x = Ezmax2_t1_x;
%     Ezmax2_t1_x = Ez(N_x, :);
% 
%     Ez1_t2_y = Ez1_t1_y;
%     Ez1_t1_y = Ez(:,1);
%     Ez2_t2_y = Ez2_t1_y;
%     Ez2_t1_y = Ez(:, 2);
% 
%     Ezmax_t2_y = Ezmax_t1_y;
%     Ezmax_t1_y = Ez(:, N_y+1);
%     Ezmax2_t2_y = Ezmax2_t1_y;
%     Ezmax2_t1_y = Ez(:, N_y);
% 
% 
%     for i = 2:N_x
%         for j = 2:N_y
%             Ez(i, j) = Ca(i, j)*Ez(i, j) + Cb(i, j)*(Hy(i, j)-Hy(i-1, j) ...
%                 + Hx(i,j-1) - Hx(i, j));
% 
%             Ez(1, j) = -Ez2_t2_x(j)-((dx-c0*dt)/(dx+c0*dt))*(Ez(2, j) + ...
%                 Ez1_t2_x(j)) + (2*dx/(dx+c0*dt))*(Ez1_t1_x(j)+Ez2_t1_x(j)) + ...
%                 (dx*(c0*dt)^2/(2*dy^2*(dx+c0*dt)))*(Ez1_t1_x(j+1)-2*Ez1_t1_x(j) + ...
%                 Ez1_t1_x(j-1)+ Ez2_t1_x(j+1)-2*Ez2_t1_x(j) +Ez2_t1_x(j-1));
% 
%             Ez(N_x+1, j) = -Ezmax2_t2_x(j)-((dx-c0*dt)/(dx+c0*dt))*(Ez(N_x, j) + ...
%                 Ezmax_t2_x(j)) + (2*dx/(dx+c0*dt))*(Ezmax_t1_x(j)+Ezmax2_t1_x(j)) + ...
%                 (dx*(c0*dt)^2/(2*dy^2*(dx+c0*dt)))*(Ezmax_t1_x(j+1)-2*Ezmax_t1_x(j) + ...
%                 Ezmax_t1_x(j-1)+ Ezmax2_t1_x(j+1)-2*Ezmax2_t1_x(j) +Ezmax2_t1_x(j-1));  
%         end
% 
%         Ez(i, 1) = -Ez2_t2_y(i)-((dx-c0*dt)/(dx+c0*dt))*(Ez(i, 2) + ...
%             Ez1_t2_y(i)) + (2*dx/(dx+c0*dt))*(Ez1_t1_y(i)+Ez2_t1_y(i)) + ...
%             (dx*(c0*dt)^2/(2*dy^2*(dx+c0*dt)))*(Ez1_t1_y(i+1)-2*Ez1_t1_y(i) + ...
%             Ez1_t1_y(i-1)+ Ez2_t1_y(i+1)-2*Ez2_t1_y(i) +Ez2_t1_y(i-1));
% 
%         Ez(i, N_y+1) = -Ezmax2_t2_y(i)-((dx-c0*dt)/(dx+c0*dt))*(Ez(i, N_y) + ...
%             Ezmax_t2_y(i)) + (2*dx/(dx+c0*dt))*(Ezmax_t1_y(i)+Ezmax2_t1_y(i)) + ...
%             (dx*(c0*dt)^2/(2*dy^2*(dx+c0*dt)))*(Ezmax_t1_y(i+1)-2*Ezmax_t1_y(i) + ...
%             Ezmax_t1_y(i-1)+ Ezmax2_t1_y(i+1)-2*Ezmax2_t1_y(i) +Ezmax2_t1_y(i-1));
%     end
%     Ez(1, 1) = Ez2_t1_x(2) - (dx-c0*dt)/(dx+c0*dt)*(Ez(2, 1) - Ez1_t1_x(1));
%     Ez(1, N_y+1) = Ezmax2_t1_x(2) - (dx-c0*dt)/(dx+c0*dt)*(Ez(2, N_y+1) - Ezmax_t1_x(1));
%     Ez(N_x+1, 1) = Ez2_t1_y(N_x) - (dx-c0*dt)/(dx+c0*dt)*(Ez(N_x+1, 2) - Ez1_t1_y(N_x+1));
%     Ez(N_x+1, N_y+1) = Ezmax2_t1_y(N_x) - (dx-c0*dt)/(dx+c0*dt)*(Ez(N_x+1, N_y) - Ezmax_t1_y(N_x+1));
% 
%     Ez(N_x/2, N_y/2) = sin(2*pi*f0*t);
% 
%     for i = 1:N_x+1
%         for j = 1:N_y
%             Hx(i, j) = Hx(i, j) - Da*(Ez(i, j+1) - Ez(i, j));
%         end
%     end
% 
%     for i = 1:N_x
%         for j = 1:N_y+1
%             Hy(i, j) = Hy(i, j) + Db*(Ez(i+1, j) - Ez(i, j));
%         end
%     end
%     
%     surf(Ez);
%     colormap default;
%     drawnow;
%    
% end


% With PML
Npml = 8;
pow = 2;
R = 10^(-6);

Hx_pml_h = zeros(Npml, N_y, 2);
Hy_pml_h = zeros(Npml, N_y, 2);
Ezx_pml_h = zeros(Npml, N_y, 2);
Ezy_pml_h = zeros(Npml, N_y, 2);

Hx_pml_v = zeros(N_x, Npml, 2);
Hy_pml_v = zeros(N_x, Npml, 2);
Ezx_pml_v = zeros(N_x, Npml, 2);
Ezy_pml_v = zeros(N_x, Npml, 2);

se= -e0*c0*log(R)/(2^(pow+2)*dx*Npml^(pow+1));
for i=1:Npml
    sigmaE(i)=se*((2*i+1)^(pow+1)-(2*i-1)^(pow+1)) ;
end

sh=m0/e0*se;
for i=1:Npml
    sigmaHx(i)=sh*((2*i+1)^(pow+1)-(2*i-1)^(pow+1));
    sigmaHy(i)=sh*((2*(i+0.5)+1)^(pow+1)-(2*(i+0.5)-1)^(pow+1));
end


sigmaE=fliplr(sigmaE); 
sigmaHx=fliplr(sigmaHx);
sigmaHy=fliplr(sigmaHy);

Ca_pml = exp(1).^(-sigmaE.*dt/e0);
Cb_pml = (1-Ca_pml)./(sigmaE.*dx);
Dax_pml = exp(1).^( -sigmaHx.*dt./m0);
Day_pml = exp(1).^( -sigmaHy.*dt./m0);
Dbx_pml = (1-Dax_pml)./(sigmaHx.*dx);
Dby_pml = (1-Day_pml)./(sigmaHy.*dx);

for t = 0:dt:Tmax

    for i = 1:N_x+1
        for j = 2:N_y
            if (i ~= 1) && (i ~= N_x+1)
                Ez(i, j) = Ca(i, j)*Ez(i, j) + Cb(i, j)*(Hy(i, j)-Hy(i-1, j) ...
                    + Hx(i,j-1) - Hx(i, j));
            elseif i == 1
                Ez(1, j) = Ca(1, j)*Ez(1, j) + Cb(1, j)*(Hy(1, j)-Hy_pml_h(Npml, j, 1) ...
                    + Hx(1,j-1) - Hx(1, j));
            elseif i == N_x+1
                Ez(N_x+1, j) = Ca(N_x+1, j)*Ez(N_x+1, j) + Cb(N_x+1, j)*(Hy_pml_h(1, j, 2) ...
                    - Hy(N_x, j) + Hx(N_x+1,j-1) - Hx(N_x+1, j));
            end
        end

        if (i ~= 1) && (i ~= N_x + 1)
            Ez(i, 1) = Ca(i, 1)*Ez(i, 1) + Cb(i, 1)*(Hy(i, 1)-Hy(i-1, 1) ...
                    + Hx_pml_v(i, Npml, 1) - Hx(i, 1));
            Ez(i, N_y+1) = Ca(i, N_y+1)*Ez(i, N_y+1) + Cb(i, N_y+1)*(Hy(i, N_y+1) ...
                -Hy(i-1, N_y+1) + Hx(i, N_y) - Hx_pml_v(i, 1, 2));
        end
        
    end

    Ez(N_x/2, N_y/2) = sin(2*pi*f0*t);


    % PML (-x)
    for i = 2:Npml
        for j = 2:N_y-1
            Ezx_pml_h(i, j, 1) = Ca_pml(i)*Ezx_pml_h(i, j, 1) + Cb_pml(i)*(Hy_pml_h(i, j, 1) ...
                - Hy_pml_h(i-1, j, 1));
            Ezy_pml_h(i, j, 1) = Ca_pml(i)*Ezy_pml_h(i, j, 1) + Cb_pml(i)*(Hx_pml_h(i, j-1, 1) ...
                - Hx_pml_h(i, j, 1));
        end
    end
    

    % PML (+x)
    for i = 1:Npml-1
        k = Npml-i+1;
        for j = 2:N_y-1
            if i ~= 1
                Ezx_pml_h(i, j, 2) = Ca_pml(k)*Ezx_pml_h(i, j, 2) + Cb_pml(k)*(Hy_pml_h(i, j, 2) ...
                    - Hy_pml_h(i-1, j, 2));
            else
                Ezx_pml_h(i, j, 2) = Ca_pml(k)*Ezx_pml_h(i, j, 2) + Cb_pml(k)*(Hy_pml_h(i, j, 2) ...
                    - Hy(N_x, j));
            end
            Ezy_pml_h(i, j, 2) = Ca_pml(k)*Ezy_pml_h(i, j, 2) + Cb_pml(k)*(Hx_pml_h(i, j-1, 2) ...
                - Hx_pml_h(i, j, 2));
        end
    end
  

    % PML (-y)
    for i = 2:N_x-1
        for j = 2:Npml
            Ezx_pml_v(i, j, 1) = Ca_pml(j)*Ezx_pml_v(i, j, 1) + Cb_pml(j)*(Hy_pml_v(i, j, 1) ...
                - Hy_pml_v(i-1, j, 1));
            Ezy_pml_v(i, j, 1) = Ca_pml(j)*Ezy_pml_v(i, j, 1) + Cb_pml(j)*(Hx_pml_v(i, j-1, 1) ...
                - Hx_pml_v(i, j, 1));
        end
    end
    

    % PML (+y)
    for i = 2:N_x-1
        for j = 1:Npml-1
            k = Npml-j+1;
            Ezx_pml_v(i, j, 2) = Ca_pml(k)*Ezx_pml_v(i, j, 2) + Cb_pml(k)*(Hy_pml_v(i, j, 2) ...
                - Hy_pml_v(i-1, j, 2));
            if j == 1
                Ezy_pml_v(i, j, 2) = Ca_pml(k)*Ezy_pml_v(i, j, 2) + Cb_pml(k)*(Hx(i, N_y) ...
                    - Hx_pml_v(i, j, 2));
            else
                Ezy_pml_v(i, j, 2) = Ca_pml(k)*Ezy_pml_v(i, j, 2) + Cb_pml(k)*(Hx_pml_v(i, j-1, 2) ...
                    - Hx_pml_v(i, j, 2));
            end
        end
    end


    for i = 1:N_x+1
        for j = 1:N_y
            Hx(i, j) = Hx(i, j) - Da*(Ez(i, j+1) - Ez(i, j));
        end
    end

    for i = 1:N_x
        for j = 1:N_y+1
            Hy(i, j) = Hy(i, j) + Db*(Ez(i+1, j) - Ez(i, j));
        end
    end


    for i = 2:Npml
        for j = 2:N_y-1
            Hx_pml_h(i, j, 1) = Dax_pml(i)*Hx_pml_h(i, j, 1) + Dbx_pml(i)*(Ezx_pml_h(i, j, 1) ...
                + Ezy_pml_h(i, j, 1) - Ezx_pml_h(i, j+1, 1) - Ezy_pml_h(i, j+1, 1));
            if i ~= Npml
                Hy_pml_h(i, j, 1) = Day_pml(i)*Hy_pml_h(i, j, 1) + Dby_pml(i)*(Ezx_pml_h(i+1, j, 1) ...
                    + Ezy_pml_h(i+1, j, 1) - Ezx_pml_h(i, j, 1) - Ezy_pml_h(i, j, 1));
            else
                Hy_pml_h(i, j, 1) = Day_pml(i)*Hy_pml_h(i, j, 1) + Dby_pml(i)*(Ez(1, j) ...
                    - Ezx_pml_h(i, j, 1) - Ezy_pml_h(i, j, 1));
            end
        end
    end

    for i = 1:Npml-1
        k = Npml-i+1;
        for j = 2:N_y-1
            Hx_pml_h(i, j, 2) = Dax_pml(k)*Hx_pml_h(i, j, 2) + Dbx_pml(k)*(Ezx_pml_h(i, j, 2) ...
                + Ezy_pml_h(i, j, 2) - Ezx_pml_h(i, j+1, 2) - Ezy_pml_h(i, j+1, 2));
            Hy_pml_h(i, j, 2) = Day_pml(k)*Hy_pml_h(i, j, 2) + Dby_pml(k)*(Ezx_pml_h(i+1, j, 2) ...
                + Ezy_pml_h(i+1, j, 2) - Ezx_pml_h(i, j, 2) - Ezy_pml_h(i, j, 2));
        end
    end

    for i = 2:N_x-1
        for j = 2:Npml
            Hy_pml_v(i, j, 1) = Day_pml(j)*Hy_pml_v(i, j, 1) + Dby_pml(j)*(Ezx_pml_v(i+1, j, 1) ...
                + Ezy_pml_v(i+1, j, 1) - Ezx_pml_v(i, j, 1) - Ezy_pml_v(i, j, 1));
            if j ~= Npml
                Hx_pml_v(i, j, 1) = Dax_pml(j)*Hx_pml_v(i, j, 1) + Dbx_pml(j)*(Ezx_pml_v(i, j, 1) ...
                    + Ezy_pml_v(i, j, 1) - Ezx_pml_v(i, j+1, 1) - Ezy_pml_v(i, j+1, 1));
            else
                Hx_pml_v(i, j, 1) = Dax_pml(j)*Hx_pml_v(i, j, 1) + Dbx_pml(j)*(Ezx_pml_v(i, j, 1) ...
                    + Ezy_pml_v(i, j, 1) - Ez(i, 1));
            end
        end
    end

    for i = 2:N_x-1
        for j = 1:Npml-1
            k = Npml-j+1;
            Hx_pml_v(i, j, 2) = Dax_pml(k)*Hx_pml_v(i, j, 2) + Dbx_pml(k)*(Ezx_pml_v(i, j, 2) ...
                + Ezy_pml_v(i, j, 2) - Ezx_pml_v(i, j+1, 2) - Ezy_pml_v(i, j+1, 2));
            Hy_pml_v(i, j, 2) = Day_pml(k)*Hy_pml_v(i, j, 2) + Dby_pml(k)*(Ezx_pml_v(i+1, j, 2) ...
                + Ezy_pml_v(i+1, j, 2) - Ezx_pml_v(i, j, 2) - Ezy_pml_v(i, j, 2));
        end
    end

    
    
    surf(Ez);
    colormap default;
    drawnow;

end

