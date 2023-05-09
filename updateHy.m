function Hy = updateHy(Hy_prev, Ez, Db, N_x, N_y)
    Hy = zeros(N_x, N_y+1);
    for i = 1:N_x
        for j = 2:N_y
            Hy(i, j) = Hy_prev(i, j) + Db*(Ez(i+1, j) - Ez(i, j));
        end
    end
end