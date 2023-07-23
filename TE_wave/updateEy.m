function Ey = updateEy(Ey_prev, Hz, Ca, Cb, N_x, N_y, Npml, has_graph, gr_pos, Cb_gr, Jy)
    Ey = Ey_prev;
    if ~has_graph
        for i = Npml+1:Npml+N_x
            for j = 1:2*Npml+N_y
                Ey(i, j) = Ca(i, j)*Ey(i, j) + Cb(i, j)*(Hz(i, j)-Hz(i+1, j));
            end        
        end
    else
        for i = Npml+1:Npml+N_x
            if i ~= gr_pos
                for j = 1:2*Npml+N_y
                    Ey(i, j) = Ca(i, j)*Ey(i, j) + Cb(i, j)*(Hz(i, j)-Hz(i+1, j));
                end   
            else
                for j = 1:2*Npml+N_y
                    Ey(i, j) = Ca(i, j)*Ey(i, j) + Cb_gr*(Hz(i, j)-Hz(i+1, j)-Jy(i, j));
                end 
            end
        end
    end
end