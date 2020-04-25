function F4 = f4(t,s,e,ya,ys,r,eta,rho_i,gamma,omega,Lim,Lsm_temp,Lmi,Lms_temp,N,h11)

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
    Lms = h11*Lms_temp;
else
    Lsm = Lsm_temp;
    Lms = Lms_temp;
end

Ni = s + e + ya + ys + r;

F4 = (1 - rho_i)*gamma*e - eta*ys - omega*ys - (Lim + Lsm)*(ys/N) + (Lmi + Lms)*(1 - Ni)*(ys/N);