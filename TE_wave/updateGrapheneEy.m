function Ey = updateGrapheneEy(Ey_prev, Hz, Ca, Cb, N_x, N_y, Npml, gr_pos, ...
    Cb_gr, Jy, hasPMLy)
    Ey = Ey_prev;
    if ~hasPMLy
        for i = Npml+1:Npml+N_x
            if i ~= gr_pos
                for j = 1:2*Npml+N_y
                    Ey(i, j) = Ca(i, j)*Ey(i, j) + Cb(i, j)*(Hz(i, j)-Hz(i+1, j));
                end   
            else
                for j = 1:2*Npml+N_y
                    Ey(i, j) = Ey(i, j) + Cb_gr*(Hz(i, j)-Hz(i+1, j) - Jy(i, j));
                end 
            end
        end
    else
        for i = Npml+1:Npml+N_x
            if i ~= gr_pos
                for j = Npml+1:Npml+N_y
                    Ey(i, j) = Ca(i, j)*Ey(i, j) + Cb(i, j)*(Hz(i, j)-Hz(i+1, j));
                end   
            else
                for j = Npml+1:Npml+N_y
                    Ey(i, j) = Ey(i, j) + Cb_gr*(Hz(i, j)-Hz(i+1, j) - Jy(i, j));
                end 
            end
        end
    end
end