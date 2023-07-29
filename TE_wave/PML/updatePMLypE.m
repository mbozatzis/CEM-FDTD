function [Ex, Ey] = updatePMLypE(Ex_prev, Ey_prev, Hzx_pml, Hzy_pml, Hz, Cax_pml, ...
    Cay_pml, Cbx_pml, Cby_pml, Npml, Nx, Ny)

    Ex = Ex_prev;
    Ey = Ey_prev;
    for i = Npml+1:Ny+Npml
        for j = Npml+Nx+1:Nx+2*Npml-1
            if j ~= Npml+Nx+1
                Ex(i, j) = Cax_pml(i, j)*Ex(i, j) + Cbx_pml(i, j)*(Hzx_pml(i, j) ...
                    + Hzy_pml(i, j) - Hzx_pml(i, j-1) - Hzy_pml(i, j-1));
            else
                Ex(i, j) = Cay_pml(i, j)*Ex(i, j) + Cby_pml(i, j)*(Hzx_pml(i, j) ...
                + Hzy_pml(i, j) - Hz(i, Npml+Ny));
            end
            Ey(i, j) = Cay_pml(i, j)*Ey(i, j) + Cby_pml(i, j)*(Hzx_pml(i, j) ...
                + Hzy_pml(i, j) - Hzx_pml(i-1, j) - Hzy_pml(i-1, j));
        end
    end
end