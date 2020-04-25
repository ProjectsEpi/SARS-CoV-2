function F4 = f44(t,e,ys,N,eta,rho_i,gamma,Lim,Lsm_temp,h11)

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
else
    Lsm = Lsm_temp;
end

F4 = (1 - rho_i)*gamma*e  - eta*ys - (Lim + Lsm)*(ys/N);