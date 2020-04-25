function F6 = f6(t,s,e,ya,ys,r,Lim,Lsm_temp,Lmi,Lms_temp,mu,N,h11)

if (43 <= t) && (t <= 51)
    Lsm = h11*Lsm_temp;
    Lms = h11*Lms_temp;
else
    Lsm = Lsm_temp;
    Lms = Lms_temp;
end

Ni = s + e + ya + ys + r;

F6 = Lim + Lsm - (Lmi + Lms)*(1 - Ni) - mu*N;