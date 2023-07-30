function Hz_inc = IncFromAuxH(H_aux, Hz_inc_prev, fi, tf_region_start, tf_region_end)

    Hz_inc = Hz_inc_prev;

    % Calculate r components for each angle
    if (fi >= 0) && (fi <= pi/2)
        base_i = tf_region_start(1);
        base_j = tf_region_start(2);
    elseif (fi > pi/2) && (fi <= pi)
        base_i = tf_region_end(1);
        base_j = tf_region_start(2);
    elseif (fi > pi) && (fi <= 3*pi/2)
        base_i = tf_region_end(1);
        base_j = tf_region_end(2);
    elseif (fi > 3*pi/2) && (fi <= 2*pi)
        base_i = tf_region_start(1);
        base_j = tf_region_end(2);
    end

    % Calculate incident field Hz from Auxilary field H - xn, xp
    for j = tf_region_start(2):tf_region_end(2)
        d1 = (tf_region_start(1) - base_i)*cos(fi) + (j - base_j)*sin(fi);
        d2 = (tf_region_end(1) - base_i)*cos(fi) + (j - base_j)*sin(fi);
        d1t = d1 - ceil(d1);
        d2t = d2 - ceil(d2);

        Hz_inc(tf_region_start(1), j) = (1-d1t)*H_aux(2+ceil(d1)) + d1t*H_aux(2+ceil(d1)+1);
        Hz_inc(tf_region_end(1), j) = (1-d2t)*H_aux(2+ceil(d2)) + d2t*H_aux(2+ceil(d2)+1);
    end

    % Calculate incident field Hz from Auxilary field H - yn, yp
    for i = tf_region_start(1):tf_region_end(1)
        d1 = (i - base_i)*cos(fi) + (tf_region_start(2) - base_j)*sin(fi);
        d2 = (i - base_i)*cos(fi) + (tf_region_end(2) - base_j)*sin(fi);
        d1t = d1 - ceil(d1);
        d2t = d2 - ceil(d2);

        Hz_inc(i, tf_region_start(2)) = (1-d1t)*H_aux(2+ceil(d1)) + d1t*H_aux(2+ceil(d1)+1);
        Hz_inc(i, tf_region_end(2)) = (1-d2t)*H_aux(2+ceil(d2)) + d2t*H_aux(2+ceil(d2)+1);
    end

end