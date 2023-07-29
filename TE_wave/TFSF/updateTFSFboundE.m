function [Ex, Ey] = updateTFSFboundE(Ex_prev, Ey_prev, Hz_inc, tf_region_start, ...
    tf_region_end, dx, dy, dt) 

    Ex = Ex_prev;
    Ey = Ey_prev;
    e0 = 8.85418781762e-12;

    % Update Ey TFSF boundaries at xn
    Ey(tf_region_start(1), tf_region_start(2):tf_region_end(2)) = ...
            Ey(tf_region_start(1), tf_region_start(2):tf_region_end(2)) + ...
            dt/(e0*dy)*Hz_inc(tf_region_start(1), tf_region_start(2):tf_region_end(2));

    % Update Ey TFSF boundaries at xp
    Ey(tf_region_end(1), tf_region_start(2):tf_region_end(2)) = ...
            Ey(tf_region_end(1), tf_region_start(2):tf_region_end(2)) - ...
            dt/(e0*dy)*Hz_inc(tf_region_end(1), tf_region_start(2):tf_region_end(2));

    % Update Ex TFSF boundaries at yn

    % Update Ex TFSF boundaries at yp


end