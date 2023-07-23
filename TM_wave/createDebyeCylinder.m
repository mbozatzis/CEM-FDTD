function [e, sigma, beta_p, k_p] = createDebyeCylinder(e_prev, sigma_prev, ...
    beta_p_prev, k_p_prev, e_inf, de, tau, sigma_r, x0, y0, R, n_lambda, dt)
    e = e_prev;
    sigma = sigma_prev;
    beta_p = beta_p_prev;
    k_p = k_p_prev;

    grid_center = [floor(size(e, 1)/2), floor(size(e, 2)/2)];
    e0 = 8.85418781762e-12;
    
    i0 = grid_center(1) + x0*n_lambda;
    j0 = grid_center(2) + y0*n_lambda;
    r = R*n_lambda;

    bp_i = (e0*de*dt./tau)./(1+dt./(2.*tau));
    kp_i = (1-dt./(2.*tau))./(1+dt./(2.*tau));

    for i = 1:size(e, 1)
        for j = 1:size(e, 2)
            if sqrt((i-i0)^2 + (j-j0)^2) <= r
                e(i, j) = e_inf*e0;
                sigma(i, j) = sigma_r;

                for l = 1:size(k_p, 3)
                    k_p(i, j, l) = kp_i(l);
                    beta_p(i, j, l) = bp_i(l);
                end
            end
        end
    end
end