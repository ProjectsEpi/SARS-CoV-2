function F3 = f33(t,e,ya,N,gamma,rho_i,q_i,rho_t,eta,Lim,Lsm_temp,h11)

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
else
    Lsm = Lsm_temp;
end

F3 = rho_i*gamma*e + (1 - q_i)*(1 - rho_t)*(Lim/N) - eta*ya - (Lim + Lsm)*(ya/N);