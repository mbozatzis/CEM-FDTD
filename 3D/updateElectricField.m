function [Ex, Ey, Ez] = updateElectricField(Ex_p, Ey_p, Ez_p, Hx, Hy, Hz, ...
    Ce, Npml, Nx, Ny, Nz, dx, dy, dz)

    Ex = Ex_p;
    Ey = Ey_p;
    Ez = Ez_p;

    % Ex
    for i = Npml+1:Npml+Nx
        for j = Npml+1:Npml+Ny
            for mk = Npml+1:Npml+Nz
                Ex(i, j, mk) = Ex(i, j, mk) + (Ce(i,j,mk)/dx)*(Hz(i, j, mk) - ...
                    Hz(i, j - 1, mk) - Hy(i, j, mk) + Hy(i, j, mk-1));
            end
        end
    end

    % Ey
    for i = Npml+1:Npml+Nx
        for j = Npml+1:Npml+Ny
            for mk = Npml+1:Npml+Nz
                Ey(i, j, mk) = Ey(i, j, mk) + (Ce(i,j,mk)/dy)*(Hx(i, j, mk) - ...
                    Hx(i, j, mk-1) - Hz(i, j, mk) + Hz(i-1, j, mk));
            end
        end
    end

    % Ez
    for i = Npml+1:Npml+Nx
        for j = Npml+1:Npml+Ny
            for mk = Npml+1:Npml+Nz
                Ez(i, j, mk) = Ez(i, j, mk) + (Ce(i,j,mk)/dz)*(Hy(i, j, mk) - ...
                    Hy(i-1, j, mk) - Hx(i, j, mk) + Hx(i, j-1, mk));
            end
        end
    end

end