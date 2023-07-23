function Ex = updateEx(Ex_prev, Hz, Ca, Cb, N_x, N_y, Npml)
    Ex = Ex_prev;
    for i = Npml+1:Npml+N_x
        for j = Npml+1:Npml+N_y
            Ex(i, j) = Ex(i, j) + Cb(i, j)*(Hz(i, j+1)-Hz(i, j));
        end        
    end
end