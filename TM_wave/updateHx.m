function Hx = updateHx(Hx_prev, Ez, Da, N_x, N_y, has_PML, Npml)
    Hx = Hx_prev;
    if ~has_PML
        for i = 2:N_x
            for j = 1:N_y
                Hx(i, j) = Hx_prev(i, j) - Da*(Ez(i, j+1) - Ez(i, j));
            end
        end
    else
        for i = Npml:Npml+N_x
            for j = Npml:Npml+N_y
                Hx(i, j) = Hx_prev(i, j) - Da*(Ez(i, j+1) - Ez(i, j));
            end
        end
    end
end