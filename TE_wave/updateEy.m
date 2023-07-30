function Ey = updateEy(Ey_prev, Hz, Ca, Cb, N_x, N_y, Npml, hasPMLy)
    Ey = Ey_prev;
    if ~hasPMLy
        for i = Npml+1:Npml+N_x
            for j = 1:2*Npml+N_y
                Ey(i, j) = Ca(i, j)*Ey(i, j) + Cb(i, j)*(Hz(i, j)-Hz(i+1, j));
            end        
        end
    else
        for i = Npml+1:Npml+N_x
            for j = Npml+1:Npml+N_y
                Ey(i, j) = Ca(i, j)*Ey(i, j) + Cb(i, j)*(Hz(i, j)-Hz(i+1, j));
            end        
        end
    end
end