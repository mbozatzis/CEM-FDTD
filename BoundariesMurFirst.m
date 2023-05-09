function Ez = BoundariesMurFirst(Ez_in, i, j, E2xprev, ENxprev, E2yprev, ENyprev, ...
    N_x, N_y, dx, dt, option)

    c0 = 3*10^8;

    if option == "xp"
        Ez = E2xprev + ((c0*dt-dx)/(c0*dt+dx))*(Ez_in(2, j) - Ez_in(1, j));
    elseif option == "xn"
        Ez = ENxprev + ((c0*dt-dx)/(c0*dt+dx))*(Ez_in(N_x, j) - Ez_in(N_x+1, j));
    elseif option == "yn"
        Ez = E2yprev + ((c0*dt-dx)/(c0*dt+dx))*(Ez_in(i, 2) - Ez_in(i, 1));
    elseif option == "yp"
        Ez = ENyprev + ((c0*dt-dx)/(c0*dt+dx))*(Ez_in(i, N_y) - Ez_in(i, N_y+1));
    end
%         Ez(1, j) = E2xprev + Cm1_x*(Ez_in(2, j) - Ez_in(1, j));
%         if option == "full"
%             Ez(N_x+1, j) = ENxprev + Cm1_x*(Ez_in(N_x, j) - Ez_in(N_x+1, j));
%         end
%     end
%     for i = 1:N_x
%         if option == "full"
%             Ez(i, 1) = E2yprev + Cm1_y*(Ez_in(i, 2) - Ez_in(i, 1));
%             Ez(i, N_y+1) = ENyprev + Cm1_y*(Ez_in(i, N_y) - Ez_in(i, N_y+1));
%         end
%     end
% 
%     if option == "full"
%         Ez(1, N_y+1) = E2xprev + Cm1_x*(Ez_in(2, N_y+1) - Ez_in(1, N_y+1));
%         Ez(N_x+1, 1) = ENxprev + Cm1_x*(Ez_in(N_x, 1) - Ez_in(N_x+1, 1));
%         Ez(1, 1) = E2yprev + Cm1_y*(Ez_in(1, 2) - Ez_in(1, 1));
%         Ez(N_x+1, N_y+1) = ENyprev + Cm1_y*(Ez_in(N_x+1, N_y) - Ez_in(N_x+1, N_y+1));
%     end
end

