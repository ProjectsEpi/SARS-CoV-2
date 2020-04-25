function F3 = f3(t,s,e,ya,ys,r,gamma,rho_i,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N,h11)

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
    Lms = h11*Lms_temp;
else
    Lsm = Lsm_temp;
    Lms = Lms_temp;
end

Ni = s + e + ya + ys + r;

F3 = rho_i*gamma*e - eta*ya - omega*ya - (Lim + Lsm)*(ya/N) + (Lmi + Lms)*(1 - Ni)*(ya/N);