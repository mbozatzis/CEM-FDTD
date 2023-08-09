function [Hx, Hy, Hz] = updateMagneticField(Hx_p, Hy_p, Hz_p, Ex, Ey, Ez, ...
    Cm, Npml, Nx, Ny, Nz, dx, dy, dz)

    Hx = Hx_p;
    Hy = Hy_p;
    Hz = Hz_p;

    % Hx
    for i = Npml+1:Npml+Nx
        for j = Npml+1:Npml+Ny
            for mk = Npml+1:Npml+Nz
                Hx(i, j, mk) = Hx(i, j, mk) + (Cm/dx)*(Ey(i, j, mk+1) - ...
                    Ey(i, j, mk) - Ez(i, j+1, mk) + Ez(i, j, mk));
            end
        end
    end

    % Hy
    for i = Npml+1:Npml+Nx
        for j = Npml+1:Npml+Ny
            for mk = Npml+1:Npml+Nz
                Hy(i, j, mk) = Hy(i, j, mk) + (Cm/dy)*(Ez(i+1, j, mk) - ...
                    Ez(i, j, mk) - Ex(i, j, mk+1) + Ex(i, j, mk));
            end
        end
    end

    % Hz
    for i = Npml+1:Npml+Nx
        for j = Npml+1:Npml+Ny
            for mk = Npml+1:Npml+Nz
                Hz(i, j, mk) = Hz(i, j, mk) + (Cm/dz)*(Ex(i, j+1, mk) - ...
                    Ex(i, j, mk) - Ey(i+1, j, mk) + Ey(i, j, mk));
            end
        end
    end

end