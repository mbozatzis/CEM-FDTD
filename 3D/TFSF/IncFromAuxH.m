function [Hx_inc, Hy_inc, Hz_inc] = IncFromAuxH(H_aux, Hx_inc_p, Hy_inc_p, ...
    Hz_inc_p, fi, theta, psi, tf_region_start, tf_region_end)

    Hx_inc = Hx_inc_p;
    Hy_inc = Hy_inc_p;
    Hz_inc = Hz_inc_p;

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

            Hx_inc(tf_region_start(1), j, mk) = ((1-d1t)*H_aux(2+ceil(d1)) + ...
                d1t*H_aux(2+ceil(d1)+1))*(sin(psi)*sin(fi)+cos(psi)*cos(theta)*cos(fi));
            Hx_inc(tf_region_end(1), j,  mk) = ((1-d2t)*H_aux(2+ceil(d2)) + ...
                d2t*H_aux(2+ceil(d2)+1))*(sin(psi)*sin(fi)+cos(psi)*cos(theta)*cos(fi));

            Hy_inc(tf_region_start(1), j, mk) = ((1-d1t)*H_aux(2+ceil(d1)) + ...
                d1t*H_aux(2+ceil(d1)+1))*(-sin(psi)*cos(fi)+cos(psi)*cos(theta)*sin(fi));
            Hy_inc(tf_region_end(1), j, mk) = ((1-d2t)*H_aux(2+ceil(d2)) + ...
                d2t*H_aux(2+ceil(d2)+1))*(-sin(psi)*cos(fi)+cos(psi)*cos(theta)*sin(fi));
        
            Hz_inc(tf_region_start(1), j, mk) = ((1-d1t)*H_aux(2+ceil(d1)) + ...
                d1t*H_aux(2+ceil(d1)+1))*(-cos(psi)*sin(theta));
            Hz_inc(tf_region_end(1), j, mk) = ((1-d2t)*H_aux(2+ceil(d2)) + ...
                d2t*H_aux(2+ceil(d2)+1))*(-cos(psi)*sin(theta));
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

            Hx_inc(i, tf_region_start(2), mk) = ((1-d1t)*H_aux(2+ceil(d1)) + ...
                d1t*H_aux(2+ceil(d1)+1))*(sin(psi)*sin(fi)+cos(psi)*cos(theta)*cos(fi));
            Hx_inc(i, tf_region_end(2), mk) = ((1-d2t)*H_aux(2+ceil(d2)) + ...
                d2t*H_aux(2+ceil(d2)+1))*(sin(psi)*sin(fi)+cos(psi)*cos(theta)*cos(fi));

            Hy_inc(i, tf_region_start(2), mk) = ((1-d1t)*H_aux(2+ceil(d1)) + ...
                d1t*H_aux(2+ceil(d1)+1))*(-sin(psi)*cos(fi)+cos(psi)*cos(theta)*sin(fi));
            Hy_inc(i, tf_region_start(2), mk) = ((1-d2t)*H_aux(2+ceil(d2)) + ...
                d2t*H_aux(2+ceil(d2)+1))*(-sin(psi)*cos(fi)+cos(psi)*cos(theta)*sin(fi));
    
            Hz_inc(i, tf_region_start(2), mk) = ((1-d1t)*H_aux(2+ceil(d1)) + ...
                d1t*H_aux(2+ceil(d1)+1))*(-cos(psi)*sin(theta));
            Hz_inc(i, tf_region_end(2), mk) = ((1-d2t)*H_aux(2+ceil(d2)) + ...
                d2t*H_aux(2+ceil(d2)+1))*(-cos(psi)*sin(theta));
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

            Hx_inc(i, j, tf_region_start(3)) = ((1-d1t)*H_aux(2+ceil(d1)) + ...
                d1t*H_aux(2+ceil(d1)+1))*(sin(psi)*sin(fi)+cos(psi)*cos(theta)*cos(fi));
            Hx_inc(i, j, tf_region_end(3)) = ((1-d2t)*H_aux(2+ceil(d2)) + ...
                d2t*H_aux(2+ceil(d2)+1))*(sin(psi)*sin(fi)+cos(psi)*cos(theta)*cos(fi));
    
            Hy_inc(i, j, tf_region_start(3)) = ((1-d1t)*H_aux(2+ceil(d1)) + ...
                d1t*H_aux(2+ceil(d1)+1))*(-sin(psi)*cos(fi)+cos(psi)*cos(theta)*sin(fi));
            Hy_inc(i, j, tf_region_end(3)) = ((1-d2t)*H_aux(2+ceil(d2)) + ...
                d2t*H_aux(2+ceil(d2)+1))*(-sin(psi)*cos(fi)+cos(psi)*cos(theta)*sin(fi));

            Hz_inc(i, j, tf_region_end(3)) = ((1-d1t)*H_aux(2+ceil(d1)) + ...
                d1t*H_aux(2+ceil(d1)+1))*(-cos(psi)*sin(theta));
            Hz_inc(i, j, tf_region_end(3)) = ((1-d2t)*H_aux(2+ceil(d2)) + ...
                d2t*H_aux(2+ceil(d2)+1))*(-cos(psi)*sin(theta));
        end
    end

end