function [Ex_inc, Ey_inc, Ez_inc] = IncFromAuxE(E_aux, Ex_inc_p, Ey_inc_p, ...
    Ez_inc_p, fi, theta, psi, tf_region_start, tf_region_end)

    Ex_inc = Ex_inc_p;
    Ey_inc = Ey_inc_p;
    Ez_inc = Ez_inc_p;

    % Calculate r components for each angle
    if (theta >= 0) && (theta <= pi/2)
        if (fi >= 0) && (fi <= pi/2)
            base_i = tf_region_start(1);
            base_j = tf_region_start(2);
            base_k = tf_region_start(3);
        elseif (fi > pi/2) && (fi <= pi)
            base_i = tf_region_end(1);
            base_j = tf_region_start(2);
            base_k = tf_region_start(3);
        elseif (fi > pi) && (fi <= 3*pi/2)
            base_i = tf_region_end(1);
            base_j = tf_region_end(2);
            base_k = tf_region_start(3);
        elseif (fi > 3*pi/2) && (fi <= 2*pi)
            base_i = tf_region_start(1);
            base_j = tf_region_end(2);
            base_k = tf_region_start(3);
        end
    elseif (theta >= pi/2) && (theta <= pi)
        if (fi >= 0) && (fi <= pi/2)
            base_i = tf_region_start(1);
            base_j = tf_region_start(2);
            base_k = tf_region_end(3);
        elseif (fi > pi/2) && (fi <= pi)
            base_i = tf_region_end(1);
            base_j = tf_region_start(2);
            base_k = tf_region_end(3);
        elseif (fi > pi) && (fi <= 3*pi/2)
            base_i = tf_region_end(1);
            base_j = tf_region_end(2);
            base_k = tf_region_end(3);
        elseif (fi > 3*pi/2) && (fi <= 2*pi)
            base_i = tf_region_start(1);
            base_j = tf_region_end(2);
            base_k = tf_region_end(3);
        end
    end


    % Calculate incident field Ez from Auxilary field E - xn, xp
    for j = tf_region_start(2):tf_region_end(2)
        for mk = tf_region_start(3):tf_region_end(3)
            d1 = (tf_region_start(1) - base_i)*cos(fi)*sin(theta) + ...
                (j - base_j)*sin(fi)*sin(theta) + (mk - base_k)*cos(theta);
            d2 = (tf_region_end(1) - base_i)*cos(fi)*sin(theta) + ...
                (j - base_j)*sin(fi)*sin(theta) + (mk - base_k)*cos(theta);
            d1t = d1 - ceil(d1);
            d2t = d2 - ceil(d2);

            Ex_inc(tf_region_start(1), j, mk) = ((1-d1t)*E_aux(2+ceil(d1)) + ...
                d1t*E_aux(2+ceil(d1)+1))*(cos(psi)*sin(fi)-sin(psi)*cos(theta)*cos(fi));
            Ex_inc(tf_region_end(1), j,  mk) = ((1-d2t)*E_aux(2+ceil(d2)) + ...
                d2t*E_aux(2+ceil(d2)+1))*(cos(psi)*sin(fi)-sin(psi)*cos(theta)*cos(fi));

            Ey_inc(tf_region_start(1), j, mk) = ((1-d1t)*E_aux(2+ceil(d1)) + ...
                d1t*E_aux(2+ceil(d1)+1))*(-cos(psi)*cos(fi)-sin(psi)*cos(theta)*sin(fi));
            Ey_inc(tf_region_end(1), j, mk) = ((1-d2t)*E_aux(2+ceil(d2)) + ...
                d2t*E_aux(2+ceil(d2)+1))*(-cos(psi)*cos(fi)-sin(psi)*cos(theta)*sin(fi));
        
            Ez_inc(tf_region_start(1), j, mk) = ((1-d1t)*E_aux(2+ceil(d1)) + ...
                d1t*E_aux(2+ceil(d1)+1))*sin(psi)*sin(theta);
            Ez_inc(tf_region_end(1), j, mk) = ((1-d2t)*E_aux(2+ceil(d2)) + ...
                d2t*E_aux(2+ceil(d2)+1))*sin(psi)*sin(theta);
        end
    end

    % Calculate incident field Ez from Auxilary field E - yn, yp
    for i = tf_region_start(1):tf_region_end(1)
        for mk = tf_region_start(3):tf_region_end(3)
            d1 = (i - base_i)*cos(fi)*sin(theta) + ...
                (tf_region_start(2) - base_j)*sin(fi)*sin(theta) + (mk - base_k)*cos(theta);
            d2 = (i - base_i)*cos(fi)*sin(theta) + ...
                (tf_region_end(2) - base_j)*sin(fi)*sin(theta) + (mk - base_k)*cos(theta);
            d1t = d1 - ceil(d1);
            d2t = d2 - ceil(d2);

            Ex_inc(i, tf_region_start(2), mk) = ((1-d1t)*E_aux(2+ceil(d1)) + ...
                d1t*E_aux(2+ceil(d1)+1))*(cos(psi)*sin(fi)-sin(psi)*cos(theta)*cos(fi));
            Ex_inc(i, tf_region_end(2), mk) = ((1-d2t)*E_aux(2+ceil(d2)) + ...
                d2t*E_aux(2+ceil(d2)+1))*(cos(psi)*sin(fi)-sin(psi)*cos(theta)*cos(fi));

            Ey_inc(i, tf_region_start(2), mk) = ((1-d1t)*E_aux(2+ceil(d1)) + ...
                d1t*E_aux(2+ceil(d1)+1))*(-cos(psi)*cos(fi)-sin(psi)*cos(theta)*sin(fi));
            Ey_inc(i, tf_region_start(2), mk) = ((1-d2t)*E_aux(2+ceil(d2)) + ...
                d2t*E_aux(2+ceil(d2)+1))*(-cos(psi)*cos(fi)-sin(psi)*cos(theta)*sin(fi));
    
            Ez_inc(i, tf_region_start(2), mk) = ((1-d1t)*E_aux(2+ceil(d1)) + ...
                d1t*E_aux(2+ceil(d1)+1))*sin(psi)*sin(theta);
            Ez_inc(i, tf_region_end(2), mk) = ((1-d2t)*E_aux(2+ceil(d2)) + ...
                d2t*E_aux(2+ceil(d2)+1))*sin(psi)*sin(theta);
        end
    end


    % Calculate incident field Ez from Auxilary field E - zn, zp
    for i = tf_region_start(1):tf_region_end(1)
        for j = tf_region_start(2):tf_region_end(2)
            d1 = (i - base_i)*cos(fi)*sin(theta) + (j - base_j)*sin(fi)*sin(theta) + ...
                (tf_region_start(3) - base_k)*cos(theta);
            d2 = (i - base_i)*cos(fi)*sin(theta) + (j - base_j)*sin(fi)*sin(theta) + ...
                (tf_region_end(3) - base_k)*cos(theta);
            d1t = d1 - ceil(d1);
            d2t = d2 - ceil(d2);

            Ex_inc(i, j, tf_region_start(3)) = ((1-d1t)*E_aux(2+ceil(d1)) + ...
                d1t*E_aux(2+ceil(d1)+1))*(cos(psi)*sin(fi)-sin(psi)*cos(theta)*cos(fi));
            Ex_inc(i, j, tf_region_end(3)) = ((1-d2t)*E_aux(2+ceil(d2)) + ...
                d2t*E_aux(2+ceil(d2)+1))*(cos(psi)*sin(fi)-sin(psi)*cos(theta)*cos(fi));
    
            Ey_inc(i, j, tf_region_start(3)) = ((1-d1t)*E_aux(2+ceil(d1)) + ...
                d1t*E_aux(2+ceil(d1)+1))*(-cos(psi)*cos(fi)-sin(psi)*cos(theta)*sin(fi));
            Ey_inc(i, j, tf_region_end(3)) = ((1-d2t)*E_aux(2+ceil(d2)) + ...
                d2t*E_aux(2+ceil(d2)+1))*(-cos(psi)*cos(fi)-sin(psi)*cos(theta)*sin(fi));

            Ez_inc(i, j, tf_region_end(3)) = ((1-d1t)*E_aux(2+ceil(d1)) + ...
                d1t*E_aux(2+ceil(d1)+1))*sin(psi)*sin(theta);
            Ez_inc(i, j, tf_region_end(3)) = ((1-d2t)*E_aux(2+ceil(d2)) + ...
                d2t*E_aux(2+ceil(d2)+1))*sin(psi)*sin(theta);
        end
    end

end