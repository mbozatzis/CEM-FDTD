function Ez = updateEzDebye(Ez_prev, Hx, Hy, Jp, Ca_debye, Cb_debye, k_p, ...
    N_x, N_y, has_PML, Npml)
    Ez = Ez_prev;
    Jp_sum = (1/2)*sum((1+k_p).*Jp, 3);
    if ~has_PML
        for i = 2:N_x
            for j = 2:N_y
                Ez(i, j) = Ca_debye(i, j)*Ez_prev(i, j) + Cb_debye(i, j)*(Hy(i, j)-Hy(i-1, j) ...
                    + Hx(i,j-1) - Hx(i, j) - Jp_sum(i, j));
            end
        end
    else
        for i = Npml+1:Npml+N_x+1
            for j = Npml+1:Npml+N_y+1
                Ez(i, j) = Ca_debye(i, j)*Ez_prev(i, j) + Cb_debye(i, j)*(Hy(i, j)-Hy(i-1, j) ...
                    + Hx(i,j-1) - Hx(i, j) - Jp_sum(i, j));
            end        
        end
    end
end