function Hz = updateHz(Hz_prev, Ex, Ey, Da, Db, N_x, N_y, Npml, hasPMLy)
    Hz = Hz_prev;
    if ~hasPMLy
        for i = Npml+1:Npml+N_x
            for j = 2:2*Npml+N_y
                Hz(i, j) = Hz(i, j) + Db*(Ex(i, j)-Ex(i, j-1) ...
                    + Ey(i-1,j) - Ey(i, j));
            end        
        end
    else
        for i = Npml+1:Npml+N_x
            for j = Npml+1:Npml+N_y
                Hz(i, j) = Hz(i, j) + Db*(Ex(i, j)-Ex(i, j-1) ...
                    + Ey(i-1,j) - Ey(i, j));
            end        
        end
    end
end