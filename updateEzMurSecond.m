function Ez_next = updateEzMurSecond(Ez, Ez_prev, Hx, Hy, Ca, Cb, N_x, N_y, dx, dy, dt, c0, boundary_case)
    Ez_next = Ez;
    A = ((dx-c0*dt)/(dx+c0*dt));
    B = (2*dx/(dx+c0*dt));
    C = (dx*(c0*dt)^2/(2*dy^2*(dx+c0*dt)));
    Cm1 = (dx-c0*dt)/(dx+c0*dt);
    
    for i = 2:N_x
        for j = 2:N_y
            % Update Electric Field
            Ez_next(i, j) = Ca(i, j)*Ez(i, j) + Cb(i, j)*(Hy(i, j)-Hy(i-1, j) ...
                + Hx(i,j-1) - Hx(i, j));

           % (xp, xn) Mur second order boundaries
            Ez_next(1, j) = -Ez_prev(1, j) - A*(Ez_next(2, j) + Ez_prev(2, j))...
                + B*(Ez(1, j) + Ez(2, j)) + C*(Ez(1, j+1) - 2*Ez(1, j) + ...
                Ez(1, j-1) + Ez(2, j+1) - 2*Ez(2, j) + Ez(2, j-1));

            if boundary_case == "full"
                Ez_next(N_x+1, j) = -Ez_prev(N_x+1, j) - A*(Ez_next(N_x, j) + ...
                    Ez_prev(N_x, j)) + B*(Ez(N_x+1, j) + Ez(N_x, j)) + ...
                    C*(Ez(N_x+1, j+1)-2*Ez(N_x+1, j) + Ez(N_x+1, j-1) + ...
                    Ez(N_x, j+1) - 2*Ez(N_x, j) + Ez(N_x, j-1));
            end
        end

        % (yp, yn) Mur second order boundaries
        if boundary_case == "full"
            Ez_next(i, 1) = -Ez_prev(i, 1)-A*(Ez_next(i, 2) + Ez_prev(i, 2)) ...
                + B*(Ez(i,1)+Ez(i, 2)) + C*(Ez(i+1,1)-2*Ez(i,1) + Ez(i-1,1) ...
                + Ez(i+1, 2)-2*Ez(i, 2) + Ez(i-1, 2));

            Ez_next(i, N_y+1) = -Ez_prev(i, N_y+1)-A*(Ez_next(i, N_y) + ...
                Ez_prev(i, N_y)) + B*(Ez(i, N_y+1)+Ez(i, N_y)) + ...
                C*(Ez(i+1, N_y+1)-2*Ez(i, N_y+1) + Ez(i-1, N_y+1) + ...
                Ez(i+1, N_y)-2*Ez(i, N_y) +Ez(i-1, N_y));
        end
    end

    % Mur first order boundaries at corners
    if boundary_case == "full"
        Ez_next(1, 1) = Ez(2, 2) - Cm1*(Ez_next(2, 1) - Ez(1, 1));
        Ez_next(1, N_y+1) = Ez(N_x, 2) - Cm1*(Ez_next(2, N_y+1) - Ez(N_x+1, 1));
        Ez_next(N_x+1, 1) = Ez(N_x, 2) - Cm1*(Ez_next(N_x+1, 2) - Ez(N_x+1,1));
        Ez_next(N_x+1, N_y+1) = Ez(N_x, N_y) - Cm1*(Ez_next(N_x+1, N_y) - Ez(N_x+1, N_y+1));
    end
end