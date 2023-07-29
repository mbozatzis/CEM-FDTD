function Jy = updateGrapheneADE(Jy_prev, Ey, tau_gr, A_gr, dt)

    Jy = Jy_prev;
    Jy = (2*tau_gr-dt)/(2*tau_gr+dt)*Jy + 2*A_gr*dt/(2*tau_gr+dt)*Ey;

end