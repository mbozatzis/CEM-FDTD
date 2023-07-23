function Jp = updateDebyeADE(Jp_prev, Ez, Ez_prev, k_p, beta_p, sigma, dt, ...
    N_x, N_y, has_PML, Npml)
    Jp = Jp_prev;
    if ~has_PML
        for i = 2:N_x
            for j = 2:N_y
                if sigma(i, j) ~= 0
                    Jp(i, j, :) = (1/2)*(1+k_p(i, j, :)).*Jp(i, j, :) + ...
                        (1/(2*dt))*beta_p(i, j, :).*(Ez(i, j) - Ez_prev(i, j));
                end
            end
        end
    else
        for i = Npml+1:Npml+N_x+1
            for j = Npml+1:Npml+N_y+1
                if sigma(i, j) ~= 0
                    Jp(i, j, :) = (1/2)*(1+k_p(i, j, :)).*Jp(i, j, :) + ...
                        (1/(2*dt))*beta_p(i, j, :).*(Ez(i, j) - Ez_prev(i, j));
                end
            end        
        end
    end

end