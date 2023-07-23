function Hy = updateHy(Hy_prev, Ez, Db, N_x, N_y, has_PML, Npml)
    Hy = Hy_prev;
    if ~has_PML
        for i = 1:N_x
            for j = 2:N_y
                Hy(i, j) = Hy_prev(i, j) + Db*(Ez(i+1, j) - Ez(i, j));
            end
        end
    else
        for i = Npml:Npml+N_x
            for j = Npml:Npml+N_y
                Hy(i, j) = Hy_prev(i, j) + Db*(Ez(i+1, j) - Ez(i, j));
            end
        end
    end
end