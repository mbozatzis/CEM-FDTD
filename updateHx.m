function Hx = updateHx(Hx_prev, Ez, Da, N_x, N_y)
    Hx = zeros(N_x+1, N_y);
    for i = 2:N_x
        for j = 1:N_y
            Hx(i, j) = Hx_prev(i, j) - Da*(Ez(i, j+1) - Ez(i, j));
        end
    end
end