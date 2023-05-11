function [e, sigma] = createCylinder(e_prev, sigma_prev, e_r, sigma_r, x0, y0, R, n_lambda)
    e = e_prev;
    sigma = sigma_prev;

    grid_center = [floor(size(e, 1)/2), floor(size(e, 2)/2)];
    e0 = 8.85418781762e-12;
    
    i0 = grid_center(1) + x0*n_lambda;
    j0 = grid_center(2) + y0*n_lambda;
    r = R*n_lambda;

    for i = 1:size(e, 1)
        for j = 1:size(e, 2)
            if sqrt((i-i0)^2 + (j-j0)^2) <= r
                e(i, j) = e_r*e0;
                sigma(i, j) = sigma_r;
            end
        end
    end
end