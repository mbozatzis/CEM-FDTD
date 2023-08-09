function [Hxy, Hxz, Hyx, Hyz, Hzx, Hzy, Hx, Hy, Hz] = updateMagneticPML(Hxy_p, Hxz_p, ...
    Hyx_p, Hyz_p, Hzx_p, Hzy_p, Hx_p, Hy_p, Hz_p, Ex, Ey, Ez, Dax_pml, Day_pml, Daz_pml, ...
    Dbx_pml, Dby_pml, Dbz_pml, Npml, Nx, Ny, Nz)

    Hxy = Hxy_p;
    Hxz = Hxz_p;
    Hyx = Hyx_p;
    Hyz = Hyz_p;
    Hzx = Hzx_p;
    Hzy = Hzy_p;

    Hx = Hx_p;
    Hy = Hy_p;
    Hz = Hz_p;

    % (-x)
    for i = 2:Npml
        for j = 2:Ny+2*Npml
            for mk = 2:Nz+2*Npml
                Hxy(i, j, mk) = Dax_pml(i, j, mk)*Hxy(i, j, mk) - ...
                    Dbx_pml(i, j, mk)*(Ez(i, j+1, mk) - Ez(i, j, mk));
                Hxz(i, j, mk) = Daz_pml(i, j, mk)*Hxz(i, j, mk) + ...
                    Dbx_pml(i, j, mk)*(Ey(i, j, mk+1) - Ey(i, j, mk));
                Hyx(i, j, mk) = Day_pml(i, j, mk)*Hyx(i, j, mk) + ...
                    Dby_pml(i, j, mk)*(Ez(i+1, j, mk) - Ez(i, j, mk));
                Hyz(i, j, mk) = Day_pml(i, j, mk)*Hyz(i, j, mk) - ...
                    Dby_pml(i, j, mk)*(Ex(i, j, mk+1) - Ex(i, j, mk));
                Hzx(i, j, mk) = Daz_pml(i, j, mk)*Hzx(i, j, mk) - ...
                    Dbz_pml(i, j, mk)*(Ey(i+1, j, mk) - Ey(i, j, mk));
                Hzy(i, j, mk) = Daz_pml(i, j, mk)*Hzy(i, j, mk) + ...
                    Dbz_pml(i, j, mk)*(Ex(i, j+1, mk) - Ex(i, j, mk));

                Hx(i, j, mk) = Hxy(i, j, mk) + Hxz(i, j, mk);
                Hy(i, j, mk) = Hyx(i, j, mk) + Hyz(i, j, mk);
                Hz(i, j, mk) = Hzx(i, j, mk) + Hzy(i, j, mk);
            end
        end
    end

    % (+x)
    for i = Npml+Nx+1:Nx+2*Npml-1
        for j = Npml+1:Nx+Npml
            for mk = Npml+1:Nx+Npml
                Hxy(i, j, mk) = Dax_pml(i, j, mk)*Hxy(i, j, mk) - ...
                    Dbx_pml(i, j, mk)*(Ez(i, j+1, mk) - Ez(i, j, mk));
                Hxz(i, j, mk) = Daz_pml(i, j, mk)*Hxz(i, j, mk) + ...
                    Dbx_pml(i, j, mk)*(Ey(i, j, mk+1) - Ey(i, j, mk));
                Hyx(i, j, mk) = Day_pml(i, j, mk)*Hyx(i, j, mk) + ...
                    Dby_pml(i, j, mk)*(Ez(i+1, j, mk) - Ez(i, j, mk));
                Hyz(i, j, mk) = Day_pml(i, j, mk)*Hyz(i, j, mk) - ...
                    Dby_pml(i, j, mk)*(Ex(i, j, mk+1) - Ex(i, j, mk));
                Hzx(i, j, mk) = Daz_pml(i, j, mk)*Hzx(i, j, mk) - ...
                    Dbz_pml(i, j, mk)*(Ey(i+1, j, mk) - Ey(i, j, mk));
                Hzy(i, j, mk) = Daz_pml(i, j, mk)*Hzy(i, j, mk) + ...
                    Dbz_pml(i, j, mk)*(Ex(i, j+1, mk) - Ex(i, j, mk));

                Hx(i, j, mk) = Hxy(i, j, mk) + Hxz(i, j, mk);
                Hy(i, j, mk) = Hyx(i, j, mk) + Hyz(i, j, mk);
                Hz(i, j, mk) = Hzx(i, j, mk) + Hzy(i, j, mk);
            end
        end
    end

    % (-y)
    for i = Npml+1:Nx+Npml
        for j = 2:Npml
            for mk = Npml+1:Nx+Npml
                Hxy(i, j, mk) = Dax_pml(i, j, mk)*Hxy(i, j, mk) - ...
                    Dbx_pml(i, j, mk)*(Ez(i, j+1, mk) - Ez(i, j, mk));
                Hxz(i, j, mk) = Daz_pml(i, j, mk)*Hxz(i, j, mk) + ...
                    Dbx_pml(i, j, mk)*(Ey(i, j, mk+1) - Ey(i, j, mk));
                Hyx(i, j, mk) = Day_pml(i, j, mk)*Hyx(i, j, mk) + ...
                    Dby_pml(i, j, mk)*(Ez(i+1, j, mk) - Ez(i, j, mk));
                Hyz(i, j, mk) = Day_pml(i, j, mk)*Hyz(i, j, mk) - ...
                    Dby_pml(i, j, mk)*(Ex(i, j, mk+1) - Ex(i, j, mk));
                Hzx(i, j, mk) = Daz_pml(i, j, mk)*Hzx(i, j, mk) - ...
                    Dbz_pml(i, j, mk)*(Ey(i+1, j, mk) - Ey(i, j, mk));
                Hzy(i, j, mk) = Daz_pml(i, j, mk)*Hzy(i, j, mk) + ...
                    Dbz_pml(i, j, mk)*(Ex(i, j+1, mk) - Ex(i, j, mk));

                Hx(i, j, mk) = Hxy(i, j, mk) + Hxz(i, j, mk);
                Hy(i, j, mk) = Hyx(i, j, mk) + Hyz(i, j, mk);
                Hz(i, j, mk) = Hzx(i, j, mk) + Hzy(i, j, mk);
            end
        end
    end

    % (+y)
    for i = Npml+1:Nx+Npml
        for j = Npml+Ny+1:Ny+2*Npml-1
            for mk = Npml+1:Nx+Npml
                Hxy(i, j, mk) = Dax_pml(i, j, mk)*Hxy(i, j, mk) - ...
                    Dbx_pml(i, j, mk)*(Ez(i, j+1, mk) - Ez(i, j, mk));
                Hxz(i, j, mk) = Daz_pml(i, j, mk)*Hxz(i, j, mk) + ...
                    Dbx_pml(i, j, mk)*(Ey(i, j, mk+1) - Ey(i, j, mk));
                Hyx(i, j, mk) = Day_pml(i, j, mk)*Hyx(i, j, mk) + ...
                    Dby_pml(i, j, mk)*(Ez(i+1, j, mk) - Ez(i, j, mk));
                Hyz(i, j, mk) = Day_pml(i, j, mk)*Hyz(i, j, mk) - ...
                    Dby_pml(i, j, mk)*(Ex(i, j, mk+1) - Ex(i, j, mk));
                Hzx(i, j, mk) = Daz_pml(i, j, mk)*Hzx(i, j, mk) - ...
                    Dbz_pml(i, j, mk)*(Ey(i+1, j, mk) - Ey(i, j, mk));
                Hzy(i, j, mk) = Daz_pml(i, j, mk)*Hzy(i, j, mk) + ...
                    Dbz_pml(i, j, mk)*(Ex(i, j+1, mk) - Ex(i, j, mk));

                Hx(i, j, mk) = Hxy(i, j, mk) + Hxz(i, j, mk);
                Hy(i, j, mk) = Hyx(i, j, mk) + Hyz(i, j, mk);
                Hz(i, j, mk) = Hzx(i, j, mk) + Hzy(i, j, mk);
            end
        end
    end

    % (-z)
    for i = Npml+1:Nx+Npml
        for j = Npml+1:Ny+Npml
            for mk = 2:Npml
                Hxy(i, j, mk) = Dax_pml(i, j, mk)*Hxy(i, j, mk) - ...
                    Dbx_pml(i, j, mk)*(Ez(i, j+1, mk) - Ez(i, j, mk));
                Hxz(i, j, mk) = Daz_pml(i, j, mk)*Hxz(i, j, mk) + ...
                    Dbx_pml(i, j, mk)*(Ey(i, j, mk+1) - Ey(i, j, mk));
                Hyx(i, j, mk) = Day_pml(i, j, mk)*Hyx(i, j, mk) + ...
                    Dby_pml(i, j, mk)*(Ez(i+1, j, mk) - Ez(i, j, mk));
                Hyz(i, j, mk) = Day_pml(i, j, mk)*Hyz(i, j, mk) - ...
                    Dby_pml(i, j, mk)*(Ex(i, j, mk+1) - Ex(i, j, mk));
                Hzx(i, j, mk) = Daz_pml(i, j, mk)*Hzx(i, j, mk) - ...
                    Dbz_pml(i, j, mk)*(Ey(i+1, j, mk) - Ey(i, j, mk));
                Hzy(i, j, mk) = Daz_pml(i, j, mk)*Hzy(i, j, mk) + ...
                    Dbz_pml(i, j, mk)*(Ex(i, j+1, mk) - Ex(i, j, mk));

                Hx(i, j, mk) = Hxy(i, j, mk) + Hxz(i, j, mk);
                Hy(i, j, mk) = Hyx(i, j, mk) + Hyz(i, j, mk);
                Hz(i, j, mk) = Hzx(i, j, mk) + Hzy(i, j, mk);
            end
        end
    end

    % (+z)
    for i = Npml+1:Nx+Npml
        for j = Npml+1:Ny+Npml
            for mk = Npml+Nz+1:Nz+2*Npml-1
                Hxy(i, j, mk) = Dax_pml(i, j, mk)*Hxy(i, j, mk) - ...
                    Dbx_pml(i, j, mk)*(Ez(i, j+1, mk) - Ez(i, j, mk));
                Hxz(i, j, mk) = Daz_pml(i, j, mk)*Hxz(i, j, mk) + ...
                    Dbx_pml(i, j, mk)*(Ey(i, j, mk+1) - Ey(i, j, mk));
                Hyx(i, j, mk) = Day_pml(i, j, mk)*Hyx(i, j, mk) + ...
                    Dby_pml(i, j, mk)*(Ez(i+1, j, mk) - Ez(i, j, mk));
                Hyz(i, j, mk) = Day_pml(i, j, mk)*Hyz(i, j, mk) - ...
                    Dby_pml(i, j, mk)*(Ex(i, j, mk+1) - Ex(i, j, mk));
                Hzx(i, j, mk) = Daz_pml(i, j, mk)*Hzx(i, j, mk) - ...
                    Dbz_pml(i, j, mk)*(Ey(i+1, j, mk) - Ey(i, j, mk));
                Hzy(i, j, mk) = Daz_pml(i, j, mk)*Hzy(i, j, mk) + ...
                    Dbz_pml(i, j, mk)*(Ex(i, j+1, mk) - Ex(i, j, mk));

                Hx(i, j, mk) = Hxy(i, j, mk) + Hxz(i, j, mk);
                Hy(i, j, mk) = Hyx(i, j, mk) + Hyz(i, j, mk);
                Hz(i, j, mk) = Hzx(i, j, mk) + Hzy(i, j, mk);
            end
        end
    end

end