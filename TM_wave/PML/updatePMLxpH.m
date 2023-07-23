function [Hx, Hy] = updatePMLxpH(Hx_prev, Hy_prev, Ez, Ezx_pml, Ezy_pml, Dax_pml, ...
    Dbx_pml, Day_pml, Dby_pml, Npml, N_x, N_y)
    Hx = Hx_prev;
    Hy = Hy_prev;
    for i = Npml+N_x+1:N_x+2*Npml-1
        for j = 2:N_y+2*Npml-1
            Hx(i, j) = Dax_pml(i, j)*Hx(i, j) + Dbx_pml(i, j)*(Ezx_pml(i, j) ...
                + Ezy_pml(i, j) - Ezx_pml(i, j+1) - Ezy_pml(i, j+1));
            Hy(i, j) = Day_pml(i, j)*Hy(i, j) + Dby_pml(i, j)*(Ezx_pml(i+1, j) ...
                + Ezy_pml(i+1, j) - Ezx_pml(i, j) - Ezy_pml(i, j));
        end
    end
end