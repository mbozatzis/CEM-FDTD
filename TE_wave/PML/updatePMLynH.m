function [Hzx_pml, Hzy_pml, Hz] = updatePMLynH(Hzx_pml_prev, Hzy_pml_prev, ...
    Hz_prev, Ex, Ey, Da_pml, Db_pml, Npml, Nx, Ny)

    Hz = Hz_prev;
    Hzx_pml = Hzx_pml_prev;
    Hzy_pml = Hzy_pml_prev;
    for i = Npml+1:Nx+Npml
        for j = 2:Npml
            Hzx_pml(i, j) = Da_pml(i, j)*Hzx_pml(i, j) + Db_pml(i, j)*(Ey(i+1, j) ...
                - Ey(i, j));
            Hzy_pml(i, j) = Da_pml(i, j)*Hzy_pml(i, j) + Db_pml(i, j)*(Ex(i, j+1) ...
                - Ex(i, j));
            Hz(i, j) = Hzx_pml(i, j) + Hzy_pml(i, j);
        end
    end 
    
end