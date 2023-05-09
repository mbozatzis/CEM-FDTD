function [eez, hhx, hhy] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options) 
    
    x0_nl = cylinder_options(1);
    R_nl = cylinder_options(2);
    e_r = cylinder_options(3);
    sigma_r = cylinder_options(4);
    y0_nl = cylinder_options(5);

    Xm_nl = simulation_options(1);
    Ym_nl = simulation_options(2);
    Tn = simulation_options(3);
    f0 = simulation_options(4);

    Npml = PML_options(1);
    pow = PML_options(2);
    Rpml = PML_options(3);


    e0=8.85418781762e-12;
    m0=1.256637061436e-06;
    c0 = 3*10^8;

    lambda0 = c0/f0;

    X_m = Xm_nl*lambda0;
    Y_m = Ym_nl*lambda0;

    % Cylinder position and characteristics
    x0 = x0_nl*lambda0;
    R = R_nl*lambda0;

    dx = lambda0/Xm_nl;
    dy = lambda0/Ym_nl;
    N_x = ceil(X_m/dx);
    N_y = ceil(Y_m/dy);

    Tmax = Tn/f0;
    dt = 0.9*(lambda0/16)/(sqrt(2)*c0);
    N_T = ceil(Tmax/dt);

    % TM: Ez, Hx, Hy
    Ez = zeros(N_x+1, N_y+1);
    Hx = zeros(N_x+1, N_y);
    Hy = zeros(N_x, N_y+1);

    e = e0*ones(N_x+1, N_y+1);
    sigma = zeros(N_x+1, N_y+1);


    % Create the cylinder
    n_lambda = lambda0/dx;
    i0 = N_x/2 - x0_nl*n_lambda;
    j0 = N_y/2 - y0_nl*n_lambda;
    r = R_nl*n_lambda;
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



    % PML Parameters    
    Hx_pml_h = zeros(Npml, N_y, 2);
    Hy_pml_h = zeros(Npml, N_y, 2);
    Ezx_pml_h = zeros(Npml, N_y, 2);
    Ezy_pml_h = zeros(Npml, N_y, 2);
    
    Hx_pml_v = zeros(N_x, Npml, 2);
    Hy_pml_v = zeros(N_x, Npml, 2);
    Ezx_pml_v = zeros(N_x, Npml, 2);
    Ezy_pml_v = zeros(N_x, Npml, 2);
    
    
    
    se= -e0*c0*log(Rpml)/(2^(pow+2)*dx*Npml^(pow+1));
    
    for i=1:Npml
        sigmaE(i)=se*((2*i+1)^(pow+1)-(2*i-1)^(pow+1));
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

    l = 1;
    if boundary == "No-boundary"
        for t = 0:dt:Tmax
            % Update Electric Field
            Ez = updateEz(Ez, Hx, Hy, Ca, Cb, N_x, N_y);

            % Update source
            Ez(N_x/2, N_y/2) = sin(2*pi*f0*t);
    
            % Update Magnetic Fields
            Hx = updateHx(Hx, Ez, Da, N_x, N_y);
            Hy = updateHy(Hy, Ez, Db, N_x, N_y);

            % Save the field values
            eez(:, :, l) = Ez;
            hhx(:, :, l) = Hx;
            hhy(:, :, l) = Hy;
            l = l + 1;
        end

    elseif boundary == "Mur-first-order"
        for t = 0:dt:Tmax         
            % Update Electrical Field with First order Mur ABC
            Ez = updateEzMurFirst(Ez, Hx, Hy, Ca, Cb, N_x, N_y, dx, dt, c0, boundary_case);

            % Update source
            Ez(N_x/2, N_y/2) = sin(2*pi*f0*t);

            % Update Magnetic Fields
            Hx = updateHx(Hx, Ez, Da, N_x, N_y);
            Hy = updateHy(Hy, Ez, Db, N_x, N_y);

            % Save the field values
            eez(:, :, l) = Ez;
            hhx(:, :, l) = Hx;
            hhy(:, :, l) = Hy;
            l = l + 1;
        end

    elseif boundary == "Mur-second-order"
        for t = 0:dt:Tmax 
            % Update Electrical Field with Second Order Mur ABC
            if l ~= 1
                Ez_prev = eez(:, :, l-1);
                Ez = updateEzMurSecond(Ez, Ez_prev, Hx, Hy, Ca, Cb, N_x, N_y, dx, dy, dt, c0, boundary_case);
            else
                Ez = updateEzMurSecond(Ez, Ez, Hx, Hy, Ca, Cb, N_x, N_y, dx, dy, dt, c0, boundary_case);
            end

            % Update source
            Ez(N_x/2, N_y/2) = sin(2*pi*f0*t);

            % Update Magnetic Fields
            Hx = updateHx(Hx, Ez, Da, N_x, N_y);
            Hy = updateHy(Hy, Ez, Db, N_x, N_y);

            % Save the field values
            eez(:, :, l) = Ez;
            hhx(:, :, l) = Hx;
            hhy(:, :, l) = Hy;
            l = l + 1;
        end

    elseif boundary == "PML"
        for t = 0:dt:Tmax
            
            % Save previous Electric field values
            E2yprev = Ez(1, 2);
            ENyprev = Ez(N_x+1, N_y);
            E2xprev = Ez(2, N_y+1);
            ENxprev = Ez(N_x, 1);

            % Update electric field with PML
            for i = 1:N_x+1
                for j = 2:N_y
                    if (i ~= 1) && (i ~= N_x+1)
                        Ez(i, j) = Ca(i, j)*Ez(i, j) + Cb(i, j)*(Hy(i, j)-Hy(i-1, j) ...
                            + Hx(i,j-1) - Hx(i, j));
                    elseif i == 1
                        Ez(1, j) = Ca(1, j)*Ez(1, j) + Cb(1, j)*(Hy(1, j)-Hy_pml_h(Npml, j, 1) ...
                            + Hx(1,j-1) - Hx(1, j));
                    elseif (i == N_x+1) && (boundary_case == "full")
                        Ez(N_x+1, j) = Ca(N_x+1, j)*Ez(N_x+1, j) + Cb(N_x+1, j)*(Hy_pml_h(1, j, 2) ...
                            - Hy(N_x, j) + Hx(N_x+1,j-1) - Hx(N_x+1, j));
                    end
                end
        
                if (i ~= 1) && (i ~= N_x + 1) && (boundary_case == "full")
                    Ez(i, 1) = Ca(i, 1)*Ez(i, 1) + Cb(i, 1)*(Hy(i, 1)-Hy(i-1, 1) ...
                            + Hx_pml_v(i, Npml, 1) - Hx(i, 1));
                    Ez(i, N_y+1) = Ca(i, N_y+1)*Ez(i, N_y+1) + Cb(i, N_y+1)*(Hy(i, N_y+1) ...
                        -Hy(i-1, N_y+1) + Hx(i, N_y) - Hx_pml_v(i, 1, 2));
                end
                
            end
        
            % Mur first order boundaries at corners
            if boundary_case == "full"
                Ez(1, N_y+1) = E2xprev + ((c0*dt-dx)/(c0*dt+dx))*(Ez(2, N_y+1) - Ez(1, N_y+1));
                Ez(N_x+1, 1) = ENxprev + ((c0*dt-dx)/(c0*dt+dx))*(Ez(N_x, 1) - Ez(N_x+1, 1));
                Ez(1, 1) = E2yprev + ((c0*dt-dx)/(c0*dt+dx))*(Ez(1, 2) - Ez(1, 1));
                Ez(N_x+1, N_y+1) = ENyprev + ((c0*dt-dx)/(c0*dt+dx))*(Ez(N_x+1, N_y) - Ez(N_x+1, N_y+1));
            end
        
            % Update source
            Ez(N_x/2, N_y/2) = sin(2*pi*f0*t);       
        
            % Update PML (-x) Electric
            for i = 2:Npml
                for j = 2:N_y-1
                    k = i;
                    Ezx_pml_h(i, j, 1) = Ca_pml(k)*Ezx_pml_h(i, j, 1) + Cb_pml(k)*(Hy_pml_h(i, j, 1) ...
                        - Hy_pml_h(i-1, j, 1));
                    Ezy_pml_h(i, j, 1) = Ca_pml(k)*Ezy_pml_h(i, j, 1) + Cb_pml(k)*(Hx_pml_h(i, j-1, 1) ...
                        - Hx_pml_h(i, j, 1));
                end
            end
            
        
            % Update PML (+x) Electric
            for i = 1:Npml-1
                k = Npml-i+1;
                for j = 2:N_y-1
                    if i ~= 1
                        Ezx_pml_h(i, j, 2) = Ca_pml(k)*Ezx_pml_h(i, j, 2) + Cb_pml(k)*(Hy_pml_h(i, j, 2) ...
                            - Hy_pml_h(i-1, j, 2));
                    else
                        Ezx_pml_h(1, j, 2) = Ca_pml(Npml)*Ezx_pml_h(1, j, 2) + Cb_pml(Npml)*(Hy_pml_h(1, j, 2) ...
                            - Hy(N_x, j));
                    end
                    Ezy_pml_h(i, j, 2) = Ca_pml(k)*Ezy_pml_h(i, j, 2) + Cb_pml(k)*(Hx_pml_h(i, j-1, 2) ...
                        - Hx_pml_h(i, j, 2));
                end
            end
          
        
            % Update PML (-y) Electric
            for i = 2:N_x-1
                for j = 2:Npml
                    k = j;
                    Ezx_pml_v(i, j, 1) = Ca_pml(k)*Ezx_pml_v(i, j, 1) + Cb_pml(k)*(Hy_pml_v(i, j, 1) ...
                        - Hy_pml_v(i-1, j, 1));
                    Ezy_pml_v(i, j, 1) = Ca_pml(k)*Ezy_pml_v(i, j, 1) + Cb_pml(k)*(Hx_pml_v(i, j-1, 1) ...
                        - Hx_pml_v(i, j, 1));
                end
            end
            
        
            % Update PML (+y) Electric
            for i = 2:N_x-1
                for j = 2:Npml
                    k = Npml-j+1;
                    Ezx_pml_v(i, j, 2) = Ca_pml(k)*Ezx_pml_v(i, j, 2) + Cb_pml(k)*(Hy_pml_v(i, j, 2) ...
                        - Hy_pml_v(i-1, j, 2));
                    if j ~= Npml
                        Ezy_pml_v(i, j, 2) = Ca_pml(k)*Ezy_pml_v(i, j, 2) + Cb_pml(k)*(Hx_pml_v(i, j-1, 2) ...
                            - Hx_pml_v(i, j, 2));
                    else
                        Ezy_pml_v(i, 1, 2) = Ca_pml(Npml)*Ezy_pml_v(i, 1, 2) + Cb_pml(Npml)*(Hx(i, N_y) ...
                            - Hx_pml_v(i, 1, 2));
                    end
                end
            end
        
            % Update Magnetic Field
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
        
            
            % Update PML (-x) Magnetic
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
        
            % Update PML (+x) Magnetic
            for i = 1:Npml-1
                k = Npml-i+1;
                for j = 2:N_y-1
                    Hx_pml_h(i, j, 2) = Dax_pml(k)*Hx_pml_h(i, j, 2) + Dbx_pml(k)*(Ezx_pml_h(i, j, 2) ...
                        + Ezy_pml_h(i, j, 2) - Ezx_pml_h(i, j+1, 2) - Ezy_pml_h(i, j+1, 2));
                    Hy_pml_h(i, j, 2) = Day_pml(k)*Hy_pml_h(i, j, 2) + Dby_pml(k)*(Ezx_pml_h(i+1, j, 2) ...
                        + Ezy_pml_h(i+1, j, 2) - Ezx_pml_h(i, j, 2) - Ezy_pml_h(i, j, 2));
                end
            end
        
            % Update PML (-y) Magnetic
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
        
            % Update PML (+y) Magnetic
            for i = 2:N_x-1
                for j = 1:Npml-1
                    k = Npml-j+1;
                    Hx_pml_v(i, j, 2) = Dax_pml(k)*Hx_pml_v(i, j, 2) + Dbx_pml(k)*(Ezx_pml_v(i, j, 2) ...
                        + Ezy_pml_v(i, j, 2) - Ezx_pml_v(i, j+1, 2) - Ezy_pml_v(i, j+1, 2));
                    Hy_pml_v(i, j, 2) = Day_pml(k)*Hy_pml_v(i, j, 2) + Dby_pml(k)*(Ezx_pml_v(i+1, j, 2) ...
                        + Ezy_pml_v(i+1, j, 2) - Ezx_pml_v(i, j, 2) - Ezy_pml_v(i, j, 2));
                end
            end

            eez(:, :, l) = Ez;
            hhx(:, :, l) = Hx;
            hhy(:, :, l) = Hy;
            l = l + 1;

        end
    end

        
end



