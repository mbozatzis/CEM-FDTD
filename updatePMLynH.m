function [Hx, Hy] = updatePMLynH(Hx_prev, Hy_prev, Ez, Ezx_pml, Ezy_pml, Dax_pml, ...
    Dbx_pml, Day_pml, Dby_pml, Npml, N_x, N_y)
    Hx = Hx_prev;
    Hy = Hy_prev;
    for i = Npml+1:N_x+Npml
        for j = 2:Npml
            Hy(i, j) = Day_pml(i, j)*Hy(i, j) + Dby_pml(i, j)*(Ezx_pml(i+1, j) ...
                + Ezy_pml(i+1, j) - Ezx_pml(i, j) - Ezy_pml(i, j));
            if j ~= Npml
                Hx(i, j) = Dax_pml(i, j)*Hx(i, j) + Dbx_pml(i, j)*(Ezx_pml(i, j) ...
                    + Ezy_pml(i, j) - Ezx_pml(i, j+1) - Ezy_pml(i, j+1));
            else
                Hx(i, j) = Dax_pml(i, j)*Hx(i, j) + Dbx_pml(i, j)*(Ezx_pml(i, j) ...
                    + Ezy_pml(i, j) - Ez(i, j));
            end
        end
    end
end