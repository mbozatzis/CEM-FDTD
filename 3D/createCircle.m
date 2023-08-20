function e = createCircle(e_prev, e_r, i0, j0, mk0, r)
    e = e_prev;

    e0 = 8.85418781762e-12;

    for i = 1:size(e, 1)
        for j = 1:size(e, 2)
            if sqrt((i-i0)^2 + (j-j0)^2) <= r
                e(i, j, mk0) = e_r*e0;
            end
        end
    end
end