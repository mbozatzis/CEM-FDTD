function Ez = updateEz(Ez_prev, Hx, Hy, Ca, Cb, N_x, N_y, has_PML, Npml)
    Ez = Ez_prev;
    if ~has_PML
        for i = 2:N_x
            for j = 2:N_y
                Ez(i, j) = Ca(i, j)*Ez_prev(i, j) + Cb(i, j)*(Hy(i, j)-Hy(i-1, j) ...
                    + Hx(i,j-1) - Hx(i, j));
            end
        end
    else
        for i = Npml+1:Npml+N_x+1
            for j = Npml+1:Npml+N_y+1
                    Ez(i, j) = Ca(i, j)*Ez_prev(i, j) + Cb(i, j)*(Hy(i, j)-Hy(i-1, j) ...
                        + Hx(i,j-1) - Hx(i, j));
            end        
        end
    end
end