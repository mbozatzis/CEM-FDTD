function Ex = updateEx(Ex_prev, Hz, Ca, Cb, N_x, N_y, Npml, hasPMLy)
    Ex = Ex_prev;
    if ~hasPMLy
        for i = Npml+1:Npml+N_x
            for j = 1:N_y+2*Npml-1
                Ex(i, j) = Ca(i, j)*Ex(i, j) + Cb(i, j)*(Hz(i, j+1)-Hz(i, j));
            end        
        end
    else
        for i = Npml+1:Npml+N_x
            for j = Npml+1:Npml+N_y
                Ex(i, j) = Ca(i, j)*Ex(i, j) + Cb(i, j)*(Hz(i, j+1)-Hz(i, j));
            end        
        end
    end
end