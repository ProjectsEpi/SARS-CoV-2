function F5 = f5(t,s,e,ya,ys,r,eta,omega,Lim,Lsm_temp,Lmi,Lms_temp,N,h11)

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
    Lms = h11*Lms_temp;
else
    Lsm = Lsm_temp;
    Lms = Lms_temp;
end

Ni = s + e + ya + ys + r;

F5 = eta*(ya + ys) - omega*r - (Lim + Lsm)*(r/N) + (Lmi + Lms)*(1 - Ni)*(r/N);