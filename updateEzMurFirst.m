function Ez_next = updateEzMurFirst(Ez, Hx, Hy, Ca, Cb, N_x, N_y, dx, dt, c0, boundary_case)
    Ez_next = Ez;
    Cm1 = ((c0*dt-dx)/(c0*dt+dx));
    for i = 1:N_x+1
        % Update Electric Field
        for j = 2:N_y
            if (i ~= 1) && (i ~= N_x + 1)
                Ez_next(i, j) = Ca(i, j)*Ez(i, j) + Cb(i, j)*(Hy(i, j)-Hy(i-1, j) ...
                    + Hx(i,j-1) - Hx(i, j));
            end
        end

        % (+/- y) Mur first order boundaries 
        if boundary_case == "full"
            Ez_next(i, 1) = Ez(i, 2) + Cm1*(Ez_next(i, 2) - Ez_next(i, 1));
            Ez_next(i, N_y+1) = Ez(i, N_y) + Cm1*(Ez_next(i, N_y) - Ez_next(i, N_y+1));
        end
    end

    % (+/- x) Mur first order boundaries 
    for j = 1:N_y+1
        Ez_next(1, j) = Ez(2, j) + Cm1*(Ez_next(2, j) - Ez_next(1, j));
        if boundary_case == "full"
            Ez_next(N_x+1, j) = Ez(N_x, j) + Cm1*(Ez_next(N_x, j) - Ez_next(N_x+1, j));
        end
    end

    % Mur first order boundaries at corners
    if boundary_case == "full"
        Ez_next(1, N_y+1) = Ez(2, N_y) + Cm1*(Ez_next(2, N_y+1) - Ez_next(1, N_y+1));
        Ez_next(N_x+1, 1) = Ez(N_x, 2) + Cm1*(Ez_next(N_x, 1) - Ez_next(N_x+1, 1));
        Ez_next(1, 1) = Ez(2, 2) + Cm1*(Ez_next(1, 2) - Ez_next(1, 1));
        Ez_next(N_x+1, N_y+1) = Ez(N_x, N_y) + Cm1*(Ez_next(N_x+1, N_y) - Ez_next(N_x+1, N_y+1));
    end
end