function [Exy, Exz, Eyx, Eyz, Ezx, Ezy, Ex, Ey, Ez] = updateElectricPML(Exy_p, Exz_p, ...
    Eyx_p, Eyz_p, Ezx_p, Ezy_p, Ex_p, Ey_p, Ez_p, Hx, Hy, Hz, Cax_pml, Cay_pml, Caz_pml, ...
    Cbx_pml, Cby_pml, Cbz_pml, Npml, Nx, Ny, Nz)

    Exy = Exy_p;
    Exz = Exz_p;
    Eyx = Eyx_p;
    Eyz = Eyz_p;
    Ezx = Ezx_p;
    Ezy = Ezy_p;

    Ex = Ex_p;
    Ey = Ey_p;
    Ez = Ez_p;

    % (-x)
    for i = 2:Npml
        for j = Npml+1:Nx+Npml
            for mk = Npml+1:Nx+Npml
                Exy(i, j, mk) = Cax_pml(i, j, mk)*Exy(i, j, mk) + ...
                    Cbx_pml(i, j, mk)*(Hz(i, j, mk) - Hz(i, j-1, mk));
                Exz(i, j, mk) = Cax_pml(i, j, mk)*Exz(i, j, mk) - ...
                    Cbx_pml(i, j, mk)*(Hy(i, j, mk) - Hy(i, j, mk-1));
                Eyx(i, j, mk) = Cay_pml(i, j, mk)*Eyx(i, j, mk) - ...
                    Cby_pml(i, j, mk)*(Hz(i, j, mk) - Hz(i-1, j, mk));
                Eyz(i, j, mk) = Cay_pml(i, j, mk)*Eyz(i, j, mk) + ...
                    Cby_pml(i, j, mk)*(Hx(i, j, mk) - Hx(i, j, mk-1));
                Ezx(i, j, mk) = Caz_pml(i, j, mk)*Ezx(i, j, mk) + ...
                    Cbz_pml(i, j, mk)*(Hy(i, j, mk) - Hy(i-1, j, mk));
                Ezy(i, j, mk) = Caz_pml(i, j, mk)*Ezy(i, j, mk) - ...
                    Cbz_pml(i, j, mk)*(Hx(i, j, mk) - Hx(i, j-1, mk));

                Ex(i, j, mk) = Exy(i, j, mk) + Exz(i, j, mk);
                Ey(i, j, mk) = Eyx(i, j, mk) + Eyz(i, j, mk);
                Ez(i, j, mk) = Ezx(i, j, mk) + Ezy(i, j, mk);
            end
        end
    end

    % (+x)
    for i = Npml+Nx+1:Nx+2*Npml-1
        for j = Npml+1:Nx+Npml
            for mk = Npml+1:Nx+Npml
                Exy(i, j, mk) = Cax_pml(i, j, mk)*Exy(i, j, mk) + ...
                    Cbx_pml(i, j, mk)*(Hz(i, j, mk) - Hz(i, j-1, mk));
                Exz(i, j, mk) = Cax_pml(i, j, mk)*Exz(i, j, mk) - ...
                    Cbx_pml(i, j, mk)*(Hy(i, j, mk) - Hy(i, j, mk-1));
                Eyx(i, j, mk) = Cay_pml(i, j, mk)*Eyx(i, j, mk) - ...
                    Cby_pml(i, j, mk)*(Hz(i, j, mk) - Hz(i-1, j, mk));
                Eyz(i, j, mk) = Cay_pml(i, j, mk)*Eyz(i, j, mk) + ...
                    Cby_pml(i, j, mk)*(Hx(i, j, mk) - Hx(i, j, mk-1));
                Ezx(i, j, mk) = Caz_pml(i, j, mk)*Ezx(i, j, mk) + ...
                    Cbz_pml(i, j, mk)*(Hy(i, j, mk) - Hy(i-1, j, mk));
                Ezy(i, j, mk) = Caz_pml(i, j, mk)*Ezy(i, j, mk) - ...
                    Cbz_pml(i, j, mk)*(Hx(i, j, mk) - Hx(i, j-1, mk));

                Ex(i, j, mk) = Exy(i, j, mk) + Exz(i, j, mk);
                Ey(i, j, mk) = Eyx(i, j, mk) + Eyz(i, j, mk);
                Ez(i, j, mk) = Ezx(i, j, mk) + Ezy(i, j, mk);
            end
        end
    end

    % (-y)
    for i = Npml+1:Nx+Npml
        for j = 2:Npml
            for mk = Npml+1:Nx+Npml
                Exy(i, j, mk) = Cax_pml(i, j, mk)*Exy(i, j, mk) + ...
                    Cbx_pml(i, j, mk)*(Hz(i, j, mk) - Hz(i, j-1, mk));
                Exz(i, j, mk) = Cax_pml(i, j, mk)*Exz(i, j, mk) - ...
                    Cbx_pml(i, j, mk)*(Hy(i, j, mk) - Hy(i, j, mk-1));
                Eyx(i, j, mk) = Cay_pml(i, j, mk)*Eyx(i, j, mk) - ...
                    Cby_pml(i, j, mk)*(Hz(i, j, mk) - Hz(i-1, j, mk));
                Eyz(i, j, mk) = Cay_pml(i, j, mk)*Eyz(i, j, mk) + ...
                    Cby_pml(i, j, mk)*(Hx(i, j, mk) - Hx(i, j, mk-1));
                Ezx(i, j, mk) = Caz_pml(i, j, mk)*Ezx(i, j, mk) + ...
                    Cbz_pml(i, j, mk)*(Hy(i, j, mk) - Hy(i-1, j, mk));
                Ezy(i, j, mk) = Caz_pml(i, j, mk)*Ezy(i, j, mk) - ...
                    Cbz_pml(i, j, mk)*(Hx(i, j, mk) - Hx(i, j-1, mk));

                Ex(i, j, mk) = Exy(i, j, mk) + Exz(i, j, mk);
                Ey(i, j, mk) = Eyx(i, j, mk) + Eyz(i, j, mk);
                Ez(i, j, mk) = Ezx(i, j, mk) + Ezy(i, j, mk);
            end
        end
    end

    % (+y)
    for i = Npml+1:Nx+Npml
        for j = Npml+Ny+1:Ny+2*Npml-1
            for mk = 2:Nz+2*Npml
                Exy(i, j, mk) = Cax_pml(i, j, mk)*Exy(i, j, mk) + ...
                    Cbx_pml(i, j, mk)*(Hz(i, j, mk) - Hz(i, j-1, mk));
                Exz(i, j, mk) = Cax_pml(i, j, mk)*Exz(i, j, mk) - ...
                    Cbx_pml(i, j, mk)*(Hy(i, j, mk) - Hy(i, j, mk-1));
                Eyx(i, j, mk) = Cay_pml(i, j, mk)*Eyx(i, j, mk) - ...
                    Cby_pml(i, j, mk)*(Hz(i, j, mk) - Hz(i-1, j, mk));
                Eyz(i, j, mk) = Cay_pml(i, j, mk)*Eyz(i, j, mk) + ...
                    Cby_pml(i, j, mk)*(Hx(i, j, mk) - Hx(i, j, mk-1));
                Ezx(i, j, mk) = Caz_pml(i, j, mk)*Ezx(i, j, mk) + ...
                    Cbz_pml(i, j, mk)*(Hy(i, j, mk) - Hy(i-1, j, mk));
                Ezy(i, j, mk) = Caz_pml(i, j, mk)*Ezy(i, j, mk) - ...
                    Cbz_pml(i, j, mk)*(Hx(i, j, mk) - Hx(i, j-1, mk));

                Ex(i, j, mk) = Exy(i, j, mk) + Exz(i, j, mk);
                Ey(i, j, mk) = Eyx(i, j, mk) + Eyz(i, j, mk);
                Ez(i, j, mk) = Ezx(i, j, mk) + Ezy(i, j, mk);
            end
        end
    end

    % (-z)
    for i = Npml+1:Nx+Npml
        for j = Npml+1:Ny+Npml
            for mk = 2:Npml
                Exy(i, j, mk) = Cax_pml(i, j, mk)*Exy(i, j, mk) + ...
                    Cbx_pml(i, j, mk)*(Hz(i, j, mk) - Hz(i, j-1, mk));
                Exz(i, j, mk) = Cax_pml(i, j, mk)*Exz(i, j, mk) - ...
                    Cbx_pml(i, j, mk)*(Hy(i, j, mk) - Hy(i, j, mk-1));
                Eyx(i, j, mk) = Cay_pml(i, j, mk)*Eyx(i, j, mk) - ...
                    Cby_pml(i, j, mk)*(Hz(i, j, mk) - Hz(i-1, j, mk));
                Eyz(i, j, mk) = Cay_pml(i, j, mk)*Eyz(i, j, mk) + ...
                    Cby_pml(i, j, mk)*(Hx(i, j, mk) - Hx(i, j, mk-1));
                Ezx(i, j, mk) = Caz_pml(i, j, mk)*Ezx(i, j, mk) + ...
                    Cbz_pml(i, j, mk)*(Hy(i, j, mk) - Hy(i-1, j, mk));
                Ezy(i, j, mk) = Caz_pml(i, j, mk)*Ezy(i, j, mk) - ...
                    Cbz_pml(i, j, mk)*(Hx(i, j, mk) - Hx(i, j-1, mk));

                Ex(i, j, mk) = Exy(i, j, mk) + Exz(i, j, mk);
                Ey(i, j, mk) = Eyx(i, j, mk) + Eyz(i, j, mk);
                Ez(i, j, mk) = Ezx(i, j, mk) + Ezy(i, j, mk);
            end
        end
    end

    % (+z)
    for i = Npml+1:Nx+Npml
        for j = Npml+1:Ny+Npml
            for mk = Npml+Nz+1:Nz+2*Npml-1
                Exy(i, j, mk) = Cax_pml(i, j, mk)*Exy(i, j, mk) + ...
                    Cbx_pml(i, j, mk)*(Hz(i, j, mk) - Hz(i, j-1, mk));
                Exz(i, j, mk) = Cax_pml(i, j, mk)*Exz(i, j, mk) - ...
                    Cbx_pml(i, j, mk)*(Hy(i, j, mk) - Hy(i, j, mk-1));
                Eyx(i, j, mk) = Cay_pml(i, j, mk)*Eyx(i, j, mk) - ...
                    Cby_pml(i, j, mk)*(Hz(i, j, mk) - Hz(i-1, j, mk));
                Eyz(i, j, mk) = Cay_pml(i, j, mk)*Eyz(i, j, mk) + ...
                    Cby_pml(i, j, mk)*(Hx(i, j, mk) - Hx(i, j, mk-1));
                Ezx(i, j, mk) = Caz_pml(i, j, mk)*Ezx(i, j, mk) + ...
                    Cbz_pml(i, j, mk)*(Hy(i, j, mk) - Hy(i-1, j, mk));
                Ezy(i, j, mk) = Caz_pml(i, j, mk)*Ezy(i, j, mk) - ...
                    Cbz_pml(i, j, mk)*(Hx(i, j, mk) - Hx(i, j-1, mk));

                Ex(i, j, mk) = Exy(i, j, mk) + Exz(i, j, mk);
                Ey(i, j, mk) = Eyx(i, j, mk) + Eyz(i, j, mk);
                Ez(i, j, mk) = Ezx(i, j, mk) + Ezy(i, j, mk);
            end
        end
    end

end